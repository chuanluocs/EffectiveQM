OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[29];
h q[48];
h q[46];
h q[4];
h q[7];
cx q[5],q[17];
cx q[7],q[19];
h q[15];
cx q[27],q[13];
cx q[41],q[8];
cx q[36],q[30];
h q[45];
h q[17];
h q[5];
h q[30];
h q[31];
cx q[34],q[27];
cx q[25],q[15];
cx q[49],q[29];
h q[30];
h q[47];
h q[4];
cx q[46],q[44];
h q[31];
cx q[4],q[52];
h q[26];
cx q[26],q[40];
h q[45];
h q[13];
h q[26];
h q[34];
cx q[48],q[1];
h q[1];
h q[15];
cx q[26],q[40];
cx q[32],q[42];
cx q[31],q[33];
cx q[43],q[51];
h q[44];
cx q[37],q[17];
h q[48];
h q[15];
h q[51];
cx q[47],q[45];
h q[46];
cx q[10],q[5];
h q[29];
cx q[23],q[17];
cx q[15],q[7];
h q[50];
h q[44];
h q[36];
h q[37];
h q[17];
h q[22];
cx q[9],q[47];
h q[29];
cx q[8],q[22];
cx q[29],q[37];
h q[22];
h q[34];
cx q[9],q[48];
cx q[23],q[17];
h q[16];
h q[3];
cx q[20],q[33];
h q[31];
h q[46];
cx q[42],q[50];
h q[22];
cx q[1],q[31];
h q[23];
h q[34];
h q[7];
h q[30];
cx q[3],q[34];
cx q[9],q[48];
cx q[9],q[47];
h q[1];
h q[23];
h q[44];
h q[5];
h q[42];
cx q[22],q[44];
h q[42];
h q[47];
h q[9];
cx q[16],q[50];
cx q[36],q[30];
cx q[46],q[44];
h q[7];
h q[1];
cx q[15],q[7];
cx q[5],q[23];
h q[34];
