OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[39];
h q[52];
cx q[33],q[47];
cx q[33],q[44];
h q[2];
h q[16];
h q[42];
h q[51];
cx q[52],q[43];
cx q[6],q[42];
cx q[42],q[32];
cx q[7],q[2];
h q[39];
h q[2];
h q[52];
h q[42];
h q[39];
cx q[39],q[29];
h q[52];
cx q[52],q[43];
h q[44];
cx q[16],q[51];
h q[52];
h q[7];
h q[47];
h q[42];
cx q[21],q[27];
cx q[33],q[47];
h q[42];
h q[0];
cx q[42],q[50];
h q[48];
cx q[13],q[42];
h q[13];
cx q[10],q[42];
cx q[36],q[19];
h q[46];
cx q[46],q[24];
h q[35];
h q[36];
cx q[48],q[0];
h q[13];
h q[1];
h q[19];
h q[24];
cx q[1],q[35];
h q[19];
cx q[19],q[11];
h q[10];
cx q[19],q[11];
h q[24];
h q[50];
cx q[19],q[11];
cx q[10],q[13];
h q[12];
cx q[6],q[8];
h q[25];
h q[41];
cx q[12],q[30];
h q[13];
cx q[46],q[24];
h q[47];
cx q[11],q[3];
cx q[41],q[52];
h q[15];
h q[30];
h q[13];
cx q[41],q[52];
h q[47];
cx q[25],q[37];
cx q[22],q[7];
h q[8];
h q[30];
h q[7];
h q[37];
h q[6];
cx q[15],q[47];
h q[7];
cx q[11],q[5];
h q[3];
h q[41];
cx q[42],q[50];
cx q[5],q[3];
h q[32];
h q[39];
h q[3];
h q[42];
h q[5];
h q[12];
h q[11];
cx q[41],q[52];
h q[39];
h q[8];
h q[52];
cx q[11],q[5];
cx q[12],q[39];
cx q[48],q[0];
cx q[8],q[32];
h q[3];
h q[39];
h q[48];
cx q[29],q[18];
cx q[39],q[21];