OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[12];
cx q[52],q[46];
h q[44];
cx q[30],q[25];
h q[49];
cx q[49],q[34];
cx q[30],q[25];
h q[21];
h q[12];
h q[21];
h q[49];
cx q[21],q[39];
h q[46];
h q[21];
h q[25];
cx q[48],q[39];
cx q[25],q[12];
cx q[52],q[44];
h q[34];
h q[52];
cx q[39],q[46];
h q[38];
h q[38];
cx q[36],q[50];
h q[45];
cx q[31],q[44];
h q[50];
h q[36];
cx q[39],q[31];
h q[36];
cx q[11],q[2];
h q[47];
h q[50];
h q[11];
h q[38];
h q[36];
h q[45];
cx q[39],q[31];
cx q[45],q[50];
cx q[17],q[47];
cx q[45],q[50];
h q[45];
cx q[2],q[38];
h q[39];
h q[44];
cx q[36],q[45];
h q[21];
h q[3];
h q[13];
cx q[38],q[3];
h q[38];
h q[15];
h q[33];
cx q[37],q[15];
cx q[36],q[45];
cx q[10],q[25];
cx q[21],q[29];
h q[9];
h q[37];
cx q[14],q[9];
h q[15];
cx q[15],q[33];
h q[21];
h q[29];
h q[13];
cx q[8],q[13];
h q[14];
h q[10];
cx q[34],q[8];
h q[9];
h q[13];
cx q[34],q[13];
h q[22];
cx q[8],q[41];
h q[14];
cx q[45],q[14];
h q[43];
h q[43];
h q[13];
h q[30];
h q[47];
cx q[35],q[10];
cx q[8],q[41];
h q[47];
h q[47];
h q[30];
h q[21];
cx q[21],q[29];
cx q[10],q[25];
cx q[22],q[30];
h q[29];
cx q[12],q[3];
h q[14];
cx q[10],q[25];
h q[12];
cx q[36],q[43];
h q[36];
h q[41];
h q[14];
cx q[47],q[14];
cx q[12],q[7];
h q[3];
h q[43];
cx q[3],q[7];
cx q[34],q[8];
h q[8];
h q[7];
h q[44];
cx q[44],q[31];
h q[43];
cx q[36],q[43];
h q[44];
cx q[7],q[34];
h q[3];
h q[31];
cx q[3],q[7];
h q[3];
h q[7];
cx q[15],q[33];
cx q[36],q[43];
h q[8];
h q[43];
h q[8];
cx q[15],q[40];
h q[21];
h q[21];
h q[33];
h q[16];
h q[40];
h q[26];
h q[13];
cx q[13],q[41];
cx q[16],q[24];
cx q[36],q[43];
cx q[1],q[21];
h q[40];
cx q[42],q[26];
h q[48];
h q[24];
h q[13];
h q[1];
cx q[40],q[33];
cx q[48],q[42];
h q[49];
h q[48];
h q[43];
cx q[49],q[26];
h q[1];
cx q[16],q[24];
cx q[43],q[50];
h q[21];
