OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
cx q[2],q[15];
h q[13];
cx q[50],q[6];
h q[38];
cx q[16],q[10];
cx q[30],q[7];
h q[26];
h q[10];
cx q[42],q[43];
h q[15];
cx q[34],q[26];
h q[25];
h q[2];
h q[31];
h q[44];
h q[42];
cx q[42],q[43];
cx q[13],q[43];
h q[43];
cx q[24],q[47];
cx q[50],q[6];
h q[44];
h q[33];
h q[30];
cx q[34],q[26];
cx q[33],q[24];
cx q[39],q[38];
h q[28];
cx q[2],q[15];
cx q[25],q[5];
cx q[28],q[31];
h q[33];
cx q[35],q[39];
h q[26];
h q[13];
h q[24];
h q[36];
cx q[28],q[35];
h q[26];
h q[28];
cx q[52],q[36];
h q[7];
h q[36];
cx q[33],q[44];
h q[47];
h q[10];
h q[33];
h q[10];
h q[36];
h q[38];
