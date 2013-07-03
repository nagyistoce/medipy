##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import multiprocessing
import select
import socket
import SocketServer
import time

try :
    import paramiko
except ImportError :
    logging.warning("medipy.io.dicom.SSHTunnelConnect not usable: cannot find paramiko package")

from Connection import Connection

class SSHTunnelConnection(Connection) :
    """ DICOM network connection through a SSH tunnel.
    
        The tunnel is as follows :
        
        (localhost, local_port) <-> (remote_host, 22) <-> (remote_host, remote_port)
    """
    
    def __init__(self, host, remote_port, calling_ae_title, called_ae_title, 
                 username, password, local_port=None) :
        
        if local_port is None :
            # Let the OS pick an available port by binding to port 0
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.bind(("", 0))
            local_port = s.getsockname()[1]
            s.close()

        Connection.__init__(self, "localhost", local_port, calling_ae_title, called_ae_title)
        
        pipe = multiprocessing.Pipe()
        
        self.process = multiprocessing.Process(
            target=SSHTunnelConnection._start_tunnel, 
            args=(host, remote_port, username, password, local_port, pipe))
        self.process.daemon = False
        self.process.start()
        
        # Wait for the child process to give some information about its success
        while(not pipe[0].poll(0.005)) :
            pass
        exception = pipe[0].recv()
        if exception :
            # Kill the child process and propagate the exception
            self._shutdown_tunnel()
            raise exception
        
    def __del__(self) :
        # Make sure the tunnel process is destroyed.
        self._shutdown_tunnel()
    
    @staticmethod
    def _start_tunnel(host, remote_port, username, password, local_port, pipe) :
        """ Start the tunnel to the remote host. This function is supposed to
            be executed in a multiprocessing.Process, so exceptions are caught
            and returned using pipe[1].send().
        """
        
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        try :
            client.connect(host, username=username, password=password)
        except Exception as e :
            pipe[1].send(e)
            return
        else :
            pipe[1].send(None)
        
        SSHTunnelConnection._forward_tunnel(
            local_port, host, remote_port, client.get_transport())
    
    @staticmethod
    def _forward_tunnel(local_port, remote_host, remote_port, transport):
        """ Start the forwarding server. This function is taken from the
            Paramiko distribution.
        """
        
        # this is a little convoluted, but lets me configure things for the Handler
        # object.  (SocketServer doesn't give Handlers any way to access the outer
        # server normally.)
        class SubHandler(Handler):
            chain_host = remote_host
            chain_port = remote_port
            ssh_transport = transport
        
        server = SocketServer.ThreadingTCPServer(("", local_port), SubHandler)
        server.daemon_threads = True
        server.allow_reuse_address = True
        
        server.serve_forever()
    
    def _shutdown_tunnel(self):
        """ Destroy the tunnel by terminating the tunnel process.
        """
        
        while(self.process.is_alive()) :
            self.process.terminate()
            time.sleep(0.005)

class Handler(SocketServer.BaseRequestHandler):
    """ Handler forwarding the data between a TCPServer and a Paramiko channel.
    """
    
    def handle(self):
        try:
            channel = self.ssh_transport.open_channel(
                "direct-tcpip", (self.chain_host, self.chain_port),
                self.request.getpeername())
        except Exception, e:
            return
        
        if channel is None:
            return
        
        done = False
        while not done :
            rlist = select.select([self.request, channel], [], [])[0]
            if self.request in rlist:
                data = self.request.recv(1024)
                if len(data) == 0:
                    done = True
                channel.send(data)
            if channel in rlist :
                data = channel.recv(1024)
                if len(data) == 0:
                    done = True
                self.request.send(data)
        
        channel.close()
        self.request.close()
