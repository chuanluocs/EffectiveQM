OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[0],q[9];
h q[0];
cx q[14],q[7];
cx q[7],q[13];
h q[7];
h q[13];
cx q[14],q[7];
h q[13];
h q[0];
h q[7];
cx q[10],q[0];
h q[14];
h q[7];
cx q[12],q[13];
h q[12];
h q[6];
cx q[1],q[6];
h q[1];
cx q[1],q[5];
cx q[1],q[5];
cx q[3],q[10];
h q[6];
h q[5];
h q[3];
h q[3];
h q[5];
