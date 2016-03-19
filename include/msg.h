#ifndef MSG_HPP_
#define MSG_HPP_

// Description: Definition for all the messages used in the program

// Error messages
#define ERRO_PARSING_OPTIONS "[E] Error in parsing options.\n"
#define ERRO_UNRECOGNIZED_OPTION "[E] Unrecognized option "
#define ERRO_NETWORK_BUILD_ERROR "[E] Network creation error.\n"
#define ERRO_ARCH_NOT_SUPPORTED "[E] Architecture type not supported.\n"
#define ERRO_RESOURCE_NOT_FOUND "[E] Resource can not be found.\n"
#define ERRO_NOT_HEADER_FLIT "[E] Can not perform the required action on header flit.\n"
#define ERRO_WRAP_CLOCK "[E] Simulation too long, net clock wrapped.\n"
#define ERRO_TRAFFIC_PATTEN_NOT_SUPPORTED "[E] Traffic pattern not supported.\n"
#define ERRO_PROC_PACKET "[E] Error in processing packet.\n"
#define ERRO_BAD_HOTSPOT_PARAM "[E] Bad hotspot traffic parameters.\n"
#define ERRO_BAD_TRACE_FILE "[E] Bad trace file "
#define ERRO_PARSING_TRACE_FILE "[E] Error in parsing trace file.\n"
#define ERRO_CONFIG_FILE "[E] Error in parsing config file.\n"
#define ERRO_APCG_FILE "[E] Error in parsing APCG file"
#define ERRO_PACKET_RATE "[E] packet generating rate can not be larger than one.\n"
#define ERRO_EQUATION_FILE "[E] can not generate equation file.\n"
#define ERRO_OPEN_FILE "[E] in opening file.\n"
#define ERRO_NETWORK_SIZE "[E] Only only square mesh network is supported.\n"
#define ERRO_SYSTEM_INIT "[E] System initialization failure\n"

// Warning messages
#define WARN_BUFFER_OVERFLOW "[W] Possible buffer overflow.\n"
#define WARN_PACKET_RATE "[W] Can not apply packet rate.\n"
#define WARN_EQN_NOT_SOLVABLE "[W] The generated performance equation may not be solvable.\n"
#define WARN_OPTION_CONFIG_FILE "[W] Unrecognized option in configuration file: "
#define WARN_EMPTY_TRACE_FILE "[W] Attaching an empty trace file: "

// Informational messages
#define INFO_TIMER_FIRED "[I] Routing engine timer fired."
#define INFO_CONFIG_NODE "[I] Apply configuration on node "
#define INFO_SIM_DONE "[I] Simulation completed successfully.\n"


#endif
