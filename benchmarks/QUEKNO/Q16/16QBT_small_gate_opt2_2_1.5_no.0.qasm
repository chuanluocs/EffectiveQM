OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[15],q[1];
cx q[15],q[1];
cx q[5],q[4];
h q[4];
h q[15];
h q[14];
h q[5];
h q[1];
h q[15];
cx q[14],q[1];
h q[3];
h q[13];
h q[4];
cx q[11],q[14];
cx q[3],q[2];
h q[4];
h q[3];
h q[0];
cx q[14],q[2];
cx q[0],q[12];
h q[14];
h q[11];
cx q[4],q[13];
h q[14];
cx q[1],q[5];
h q[3];
h q[14];
cx q[1],q[5];
cx q[14],q[2];
cx q[11],q[3];
h q[4];
h q[11];
h q[12];
