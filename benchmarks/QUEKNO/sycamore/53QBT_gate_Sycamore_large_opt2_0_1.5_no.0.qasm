OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
cx q[15],q[3];
h q[10];
h q[4];
h q[45];
cx q[29],q[20];
cx q[15],q[49];
h q[29];
h q[5];
h q[6];
h q[20];
cx q[29],q[20];
h q[3];
cx q[21],q[51];
cx q[50],q[45];
cx q[42],q[5];
cx q[50],q[24];
cx q[43],q[44];
h q[45];
cx q[4],q[23];
h q[42];
cx q[21],q[51];
h q[17];
h q[39];
cx q[34],q[16];
cx q[50],q[45];
h q[20];
h q[44];
h q[51];
cx q[20],q[30];
h q[35];
h q[43];
cx q[16],q[6];
h q[42];
cx q[24],q[35];
h q[44];
h q[44];
cx q[42],q[5];
h q[43];
cx q[17],q[49];
h q[17];
h q[23];
h q[29];
h q[5];
h q[51];
cx q[10],q[51];
h q[39];
h q[49];
h q[35];
h q[23];
cx q[32],q[39];