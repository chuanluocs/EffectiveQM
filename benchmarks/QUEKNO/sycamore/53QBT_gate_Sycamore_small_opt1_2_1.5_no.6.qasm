OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[10];
h q[8];
h q[44];
cx q[31],q[8];
h q[8];
h q[52];
cx q[45],q[50];
h q[15];
cx q[10],q[45];
h q[51];
cx q[21],q[52];
h q[15];
cx q[50],q[30];
cx q[44],q[41];
cx q[30],q[46];
cx q[7],q[47];
h q[45];
h q[51];
h q[30];
h q[45];
h q[30];
h q[30];
cx q[43],q[51];
cx q[8],q[15];
h q[44];
cx q[10],q[45];
h q[51];
h q[51];
h q[21];
cx q[50],q[46];
h q[45];
h q[21];
cx q[38],q[26];
h q[26];
cx q[25],q[40];
h q[26];
h q[20];
h q[48];
h q[30];
h q[30];
h q[21];
cx q[15],q[39];
h q[4];
cx q[48],q[1];
cx q[30],q[47];
h q[45];
cx q[38],q[26];
cx q[40],q[21];
cx q[20],q[4];
h q[39];
h q[25];
h q[19];
cx q[35],q[45];
h q[39];
cx q[19],q[1];
h q[19];
h q[38];
h q[21];
cx q[48],q[19];
cx q[33],q[12];
h q[38];
cx q[10],q[33];
h q[46];
h q[10];
cx q[46],q[34];
h q[38];
cx q[36],q[8];
h q[19];
h q[8];
h q[42];
cx q[10],q[33];
h q[36];
h q[36];
h q[33];
cx q[38],q[26];
h q[19];
h q[36];
cx q[42],q[33];
cx q[21],q[14];
h q[48];
cx q[21],q[14];
