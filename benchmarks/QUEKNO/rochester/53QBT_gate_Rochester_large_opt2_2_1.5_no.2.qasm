OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[48];
h q[2];
h q[16];
h q[22];
h q[27];
cx q[23],q[49];
h q[52];
cx q[52],q[1];
h q[10];
cx q[48],q[25];
h q[19];
h q[41];
cx q[31],q[33];
cx q[22],q[45];
h q[45];
h q[35];
cx q[50],q[5];
cx q[18],q[40];
h q[1];
h q[1];
h q[19];
cx q[19],q[27];
h q[36];
cx q[40],q[36];
h q[44];
h q[23];
cx q[5],q[18];
h q[40];
h q[33];
cx q[35],q[7];
h q[27];
cx q[16],q[37];
h q[1];
h q[51];
cx q[10],q[41];
h q[13];
h q[33];
h q[2];
cx q[25],q[2];
cx q[25],q[2];
cx q[51],q[44];
cx q[13],q[23];
h q[31];
h q[13];
h q[4];
cx q[52],q[12];
h q[26];
h q[7];
cx q[46],q[26];
h q[46];
cx q[51],q[44];
h q[23];
h q[7];
cx q[34],q[38];
cx q[34],q[38];
h q[13];
h q[35];
h q[12];
h q[47];
h q[44];
h q[4];
cx q[22],q[39];
h q[39];
h q[35];
cx q[13],q[23];
h q[46];
h q[10];
cx q[12],q[4];
cx q[46],q[26];
h q[30];
cx q[34],q[38];
h q[15];
h q[12];
cx q[15],q[35];
cx q[7],q[30];
h q[7];
cx q[42],q[52];
cx q[24],q[22];
h q[31];
cx q[43],q[47];
cx q[26],q[12];
h q[13];
h q[47];
h q[26];
h q[38];
cx q[31],q[10];