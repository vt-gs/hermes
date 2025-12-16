import socket
import numpy

#COSMOS TCP Connection
HOST = '100.64.109.56'  # Listen on all available interfaces
PORT = 9090

#message_payload = b"\x48\x65\x6c\x6c\x6f\x20\x57\x6f\x72\x6c\x64\x21" #"Hello World!" 
message_21 = b"\x21\x00\x00\x00\x00\x00\x11\x52\xaf\x00\x13\xcc\x78\x25\x03\xa0\x00\xa0\x00\x67\x00\x00\x00\x94"

def write_data(message):
    FEND = numpy.uint8(0xc0)
    FESC = numpy.uint8(0xdb)
    TFEND = numpy.uint8(0xdc)
    TFESC = numpy.uint8(0xdd)

    pkt_dict = {"dest_call":"",
        "dest_ssid":0,
        "src_call":"",
        "src_ssid":1,
        "control":"03",
        "pid":"F0"}  # Add FCC Callsigns

    #DATA TO CSP+AX25

    packet = bytearray()
    #--Handle Destination Call----
    for char in pkt_dict['dest_call']: packet.append(((ord(char) << 1) & 0xFE))
    #--Handle Destination SSID----
    packet.append(((pkt_dict['dest_ssid'] << 1) | 0x60))
    #--Handle Source Call----
    for char in pkt_dict['src_call']: packet.append(((ord(char) << 1) & 0xFE))
    #--Handle Source SSID----
    packet.append(((pkt_dict['src_ssid'] << 1) | 0x61))
    #--Handle Control Byte----
    packet.append(0x03) #Contol
    #--Handle PID Byte----
    packet.append(0xF0) #PID
    #--Handle Message----
    #message = b'\xa6\xa4\x47\x00\x10\x17\x19\x96' + data #Data with CSP Header
    packet.extend(message)
    # print ("AX.25 Frame: {:s}".format(packet.hex()))

    # AX25 TO KISS

    kiss = bytearray()
    kiss.append(FEND)
    kiss.append(numpy.uint8(0))
    for i,x in enumerate(packet):
        if x == FESC:
            kiss.append(FESC)
            kiss.append(TFESC)
        elif x == FEND:
            kiss.append(FESC)
            kiss.append(TFEND)
        else:
            kiss.append(numpy.uint8(x))
    kiss.append(FEND)

    return (bytes(kiss))

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.connect((HOST, PORT))
    print(f"Connected to {HOST}:{PORT}")
    resp = write_data(message_21)
    s.sendall(resp)
    print("Sent:", resp.hex())