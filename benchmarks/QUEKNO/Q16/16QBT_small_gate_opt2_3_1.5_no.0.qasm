OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
h q[10];
h q[5];
cx q[9],q[3];
h q[5];
h q[5];
h q[2];
h q[6];
cx q[9],q[3];
cx q[2],q[5];
h q[3];
cx q[6],q[8];
h q[6];
h q[2];
cx q[9],q[14];
cx q[5],q[10];
h q[9];
h q[3];
cx q[1],q[14];
h q[10];
cx q[10],q[3];
cx q[15],q[11];
h q[1];
h q[3];
cx q[9],q[14];
cx q[10],q[3];
h q[1];
h q[9];
h q[10];
cx q[1],q[14];
h q[11];
h q[11];
cx q[3],q[11];
cx q[13],q[9];
h q[9];
h q[13];
h q[3];
h q[11];
h q[11];
h q[11];
cx q[1],q[14];
cx q[3],q[14];