OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
cx q[45],q[30];
h q[25];
h q[4];
cx q[11],q[43];
cx q[10],q[37];
cx q[24],q[17];
h q[25];
cx q[49],q[20];
h q[49];
cx q[10],q[37];
h q[25];
h q[43];
cx q[25],q[32];
h q[24];
cx q[37],q[4];
cx q[34],q[24];
h q[25];
h q[34];
h q[24];
h q[34];
h q[30];
h q[43];
cx q[28],q[24];
h q[20];
h q[10];
h q[32];
cx q[37],q[4];
h q[34];
h q[6];
h q[21];
cx q[12],q[43];
h q[12];
h q[13];
cx q[11],q[12];
cx q[23],q[47];
h q[11];
h q[13];
cx q[37],q[13];
h q[12];
cx q[35],q[3];
h q[47];
cx q[6],q[29];
h q[31];
h q[46];
h q[8];
cx q[23],q[47];
h q[13];
cx q[22],q[31];
cx q[21],q[31];
cx q[8],q[46];
h q[8];
h q[3];
cx q[0],q[9];
h q[6];
h q[35];
h q[37];
cx q[21],q[7];
h q[17];
h q[8];
h q[24];
h q[46];
h q[24];
cx q[51],q[23];
cx q[44],q[45];
h q[24];
cx q[31],q[7];
h q[8];
h q[22];
cx q[8],q[46];
cx q[22],q[10];
h q[15];
cx q[24],q[45];
h q[24];
h q[24];
cx q[9],q[17];
h q[45];
h q[42];
cx q[15],q[42];
h q[9];
cx q[8],q[46];
h q[31];
cx q[45],q[23];
h q[4];
h q[0];
cx q[15],q[20];
h q[48];
cx q[6],q[1];
h q[6];
h q[20];
h q[20];
cx q[48],q[22];
h q[48];
cx q[4],q[46];
h q[6];
h q[31];
h q[22];
cx q[45],q[51];
h q[15];
h q[0];
h q[51];
cx q[15],q[20];
h q[9];
cx q[21],q[7];
cx q[0],q[9];
cx q[22],q[31];
h q[31];
cx q[48],q[31];
h q[30];
h q[6];
cx q[15],q[20];
cx q[13],q[6];
h q[48];
h q[13];
h q[23];
h q[38];
cx q[22],q[10];
cx q[39],q[38];
cx q[22],q[10];
h q[13];
h q[48];
cx q[30],q[18];
h q[20];
h q[13];
cx q[10],q[23];
h q[18];
h q[23];