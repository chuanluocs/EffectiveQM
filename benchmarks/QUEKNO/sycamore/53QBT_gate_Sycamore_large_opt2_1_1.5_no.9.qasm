OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[7];
cx q[31],q[6];
h q[31];
h q[18];
h q[7];
h q[24];
h q[24];
h q[15];
h q[48];
cx q[2],q[45];
h q[26];
cx q[18],q[16];
h q[41];
cx q[51],q[25];
cx q[3],q[50];
h q[41];
cx q[27],q[4];
cx q[48],q[0];
h q[7];
h q[35];
cx q[27],q[35];
h q[8];
h q[7];
h q[8];
cx q[26],q[30];
cx q[51],q[34];
cx q[27],q[35];
h q[25];
h q[7];
cx q[41],q[33];
cx q[40],q[41];
cx q[25],q[3];
cx q[8],q[40];
h q[51];
h q[3];
h q[24];
cx q[39],q[3];
h q[30];
cx q[11],q[24];
h q[7];
h q[24];
h q[51];
cx q[15],q[7];
h q[2];
h q[2];
h q[2];
h q[25];
cx q[25],q[50];
cx q[41],q[33];
h q[41];
h q[50];
cx q[8],q[15];
h q[15];
cx q[37],q[47];
h q[2];
h q[50];
h q[15];
cx q[8],q[15];
h q[41];
cx q[38],q[9];
h q[1];
cx q[30],q[44];
cx q[23],q[13];
h q[41];
h q[9];
h q[28];
cx q[8],q[40];
cx q[50],q[41];
h q[13];
h q[43];
h q[23];
h q[33];
cx q[16],q[15];
h q[25];
cx q[44],q[13];
h q[1];
cx q[23],q[2];
cx q[1],q[26];
cx q[19],q[1];
cx q[30],q[44];
cx q[10],q[28];
h q[16];
h q[25];
h q[16];
cx q[43],q[29];
cx q[16],q[15];
h q[15];
cx q[1],q[26];
h q[25];
h q[47];
h q[10];
h q[38];
h q[33];
h q[50];
h q[30];