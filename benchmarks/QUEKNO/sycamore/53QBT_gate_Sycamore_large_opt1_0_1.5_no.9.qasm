OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
h q[40];
h q[30];
cx q[40],q[28];
h q[23];
h q[35];
h q[23];
cx q[35],q[50];
h q[13];
cx q[10],q[30];
h q[4];
h q[18];
cx q[33],q[48];
h q[32];
h q[50];
h q[47];
h q[17];
cx q[39],q[43];
h q[23];
h q[35];
h q[21];
cx q[28],q[49];
h q[21];
cx q[37],q[8];
h q[40];
cx q[47],q[32];
h q[33];
cx q[4],q[17];
h q[32];
h q[36];
cx q[35],q[21];
h q[8];
h q[13];
cx q[45],q[2];
cx q[21],q[13];
cx q[36],q[22];
h q[50];
h q[39];
h q[4];
h q[33];
cx q[4],q[17];
cx q[18],q[37];
cx q[17],q[23];
cx q[17],q[23];
