OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
cx q[13],q[11];
h q[15];
h q[0];
cx q[4],q[18];
cx q[8],q[0];
h q[0];
cx q[2],q[15];
h q[18];
h q[15];
h q[1];
cx q[17],q[10];
cx q[11],q[4];
h q[10];
h q[17];
h q[1];
cx q[12],q[5];
cx q[18],q[1];
h q[19];
cx q[8],q[19];
h q[17];
cx q[1],q[10];
h q[4];
h q[18];
h q[10];
h q[10];
cx q[18],q[1];
h q[1];
h q[10];
