OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
cx q[40],q[47];
h q[32];
h q[8];
h q[11];
h q[23];
cx q[6],q[52];
h q[31];
h q[42];
h q[16];
h q[23];
h q[6];
cx q[39],q[23];
h q[3];
h q[31];
h q[8];
h q[36];
cx q[8],q[11];
h q[40];
cx q[11],q[3];
cx q[31],q[0];
cx q[16],q[52];
h q[6];
cx q[18],q[36];
h q[31];
cx q[42],q[32];
cx q[16],q[52];
cx q[6],q[52];
h q[52];
h q[14];
cx q[23],q[3];
h q[19];
cx q[37],q[44];
h q[9];
h q[3];
h q[14];
cx q[14],q[13];
h q[32];
cx q[12],q[9];
h q[14];
cx q[10],q[31];
h q[3];
cx q[0],q[20];
h q[14];
h q[12];
h q[14];
cx q[8],q[39];
cx q[19],q[46];
h q[28];
h q[8];
h q[19];
h q[10];
cx q[38],q[48];
h q[20];
cx q[32],q[28];
cx q[37],q[44];
h q[19];
h q[0];
cx q[14],q[13];
h q[20];
cx q[18],q[46];
h q[38];
h q[18];
cx q[14],q[3];
cx q[13],q[9];
cx q[8],q[39];
cx q[34],q[20];
h q[48];
h q[48];
h q[34];
cx q[32],q[28];
h q[28];
h q[18];
cx q[52],q[8];
h q[19];
h q[8];
h q[14];
cx q[32],q[28];
h q[39];
h q[13];
cx q[7],q[21];
cx q[38],q[48];
h q[48];
h q[20];
h q[18];
h q[46];
cx q[18],q[19];
cx q[8],q[39];
h q[7];
