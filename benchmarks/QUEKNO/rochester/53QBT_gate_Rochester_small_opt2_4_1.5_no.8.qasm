OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[32];
h q[52];
h q[22];
h q[22];
cx q[46],q[28];
cx q[20],q[30];
cx q[22],q[35];
h q[28];
cx q[34],q[52];
cx q[25],q[45];
cx q[28],q[23];
h q[24];
h q[22];
cx q[24],q[37];
h q[37];
h q[22];
h q[37];
cx q[10],q[32];
h q[20];
cx q[10],q[32];
h q[20];
h q[34];
h q[22];
h q[11];
cx q[51],q[35];
h q[34];
h q[48];
cx q[25],q[3];
h q[25];
h q[9];
cx q[33],q[19];
h q[51];
h q[11];
h q[13];
h q[25];
cx q[9],q[13];
h q[41];
cx q[41],q[37];
cx q[11],q[21];
cx q[9],q[13];
cx q[4],q[48];
h q[37];
h q[32];
h q[34];
cx q[32],q[35];
h q[48];
h q[44];
cx q[35],q[34];
cx q[44],q[46];
h q[19];
h q[35];
cx q[46],q[37];
cx q[44],q[24];
h q[46];
h q[27];
h q[13];
h q[22];
cx q[50],q[39];
cx q[24],q[46];
h q[9];
h q[24];
h q[22];
h q[27];
cx q[22],q[27];
cx q[9],q[1];
cx q[50],q[39];
h q[13];
h q[13];
cx q[38],q[8];
h q[38];
h q[50];
h q[9];
cx q[9],q[13];
h q[8];
