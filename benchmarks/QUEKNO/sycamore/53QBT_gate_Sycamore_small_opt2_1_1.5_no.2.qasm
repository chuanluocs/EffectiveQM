OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
cx q[7],q[47];
h q[48];
h q[48];
cx q[7],q[47];
h q[9];
cx q[42],q[14];
cx q[2],q[24];
h q[8];
h q[47];
cx q[24],q[9];
h q[39];
h q[24];
h q[37];
cx q[10],q[48];
h q[9];
cx q[8],q[22];
cx q[37],q[39];
h q[37];
h q[7];
cx q[19],q[6];
h q[6];
h q[9];
h q[14];
cx q[42],q[14];
h q[47];
h q[30];
cx q[8],q[5];
h q[22];
h q[5];
h q[5];
cx q[18],q[43];
h q[23];
cx q[18],q[43];
cx q[16],q[11];
h q[5];
cx q[18],q[43];
cx q[28],q[31];
h q[37];
h q[30];
h q[16];
cx q[27],q[16];
h q[31];
cx q[37],q[39];
h q[27];
h q[8];
h q[30];
h q[51];
h q[25];
cx q[51],q[30];
h q[23];
cx q[25],q[23];
cx q[5],q[22];
h q[23];