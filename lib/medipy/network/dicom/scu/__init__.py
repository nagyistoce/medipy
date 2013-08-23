from Move import Move
from scu import Echo, Get, Find, SCU

#: Query keys for each level.
keys = {
    "patient" : {
        "patient" : [(0x0010,0x0020)], # Patient ID 
        "study" : [(0x0020,0x000d)], # Study Instance UID
        "series" : [(0x0020,0x000e)], # Series Instance UID
        "image" : [(0x0008,0x0018)], # SOP Instance UID 
    },
    "study" : {
        "study" : [(0x0020,0x000d)], # Study Instance UID
        "series" : [(0x0020,0x000e)], # Series Instance UID
        "image" : [(0x0008,0x0018)], # SOP Instance UID
    }
}
