OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[34];
h q[32];
h q[39];
cx q[20],q[21];
cx q[38],q[13];
cx q[47],q[39];
cx q[4],q[52];
h q[38];
h q[38];
h q[32];
cx q[20],q[34];
h q[47];
h q[22];
h q[4];
cx q[22],q[6];
h q[4];
h q[48];
cx q[38],q[13];
h q[4];
h q[22];
h q[38];
cx q[48],q[8];
h q[13];
cx q[5],q[20];
cx q[41],q[32];
h q[49];
h q[36];
cx q[34],q[21];
cx q[47],q[39];
cx q[32],q[36];
h q[39];
cx q[49],q[32];
h q[36];
h q[39];
cx q[25],q[26];
cx q[27],q[24];
h q[34];
h q[36];
h q[21];
h q[21];
cx q[28],q[35];
h q[39];
cx q[27],q[24];
h q[27];
h q[36];
h q[8];
h q[21];
h q[8];
cx q[48],q[8];
cx q[39],q[26];
