OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[43];
h q[2];
h q[24];
h q[11];
h q[31];
cx q[20],q[46];
cx q[18],q[47];
cx q[45],q[27];
cx q[12],q[21];
cx q[31],q[24];
h q[17];
cx q[48],q[14];
h q[29];
cx q[18],q[47];
h q[24];
cx q[17],q[38];
h q[7];
h q[27];
cx q[5],q[43];
h q[27];
cx q[46],q[2];
h q[21];
cx q[4],q[1];
cx q[51],q[11];
h q[2];
cx q[16],q[35];
cx q[9],q[18];
h q[43];
cx q[41],q[52];
h q[2];
h q[35];
h q[9];
h q[9];
cx q[8],q[29];
h q[14];
h q[14];
h q[16];
h q[5];
h q[51];
h q[45];
cx q[46],q[2];
h q[46];
h q[16];
h q[2];
h q[52];
cx q[8],q[29];
cx q[11],q[7];
h q[1];
h q[9];
cx q[11],q[7];
cx q[51],q[7];
h q[36];
h q[36];
h q[24];
h q[5];
h q[16];
cx q[18],q[47];
h q[6];
h q[5];
cx q[24],q[49];
cx q[6],q[13];
h q[23];
h q[50];
cx q[26],q[39];
cx q[39],q[48];
h q[44];
h q[16];
h q[49];
cx q[32],q[27];
h q[36];
h q[37];
h q[35];
h q[7];
cx q[44],q[8];
h q[6];
h q[7];
cx q[50],q[37];
h q[50];
h q[49];
h q[43];
h q[39];
cx q[6],q[13];
h q[6];
cx q[32],q[27];
cx q[23],q[29];
h q[32];
h q[49];
cx q[40],q[31];
cx q[17],q[51];
h q[7];
h q[48];
cx q[5],q[43];
h q[50];
cx q[16],q[36];
h q[8];
h q[16];
cx q[50],q[37];
cx q[24],q[15];
cx q[16],q[35];
cx q[39],q[48];
