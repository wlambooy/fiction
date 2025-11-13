// Benchmark "ex" written by ABC on Fri Jul  4 15:14:14 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6, new_n7, new_n8, new_n9, new_n10, new_n11, new_n12, new_n13,
    new_n14, new_n15, new_n16, new_n17, new_n18, new_n19, new_n20;
  assign new_n6 = b & c;
  assign new_n7 = ~b & ~c;
  assign new_n8 = ~new_n6 & ~new_n7;
  assign new_n9 = a & ~new_n8;
  assign new_n10 = ~b & c;
  assign new_n11 = b & ~c;
  assign new_n12 = ~new_n10 & ~new_n11;
  assign new_n13 = ~a & ~new_n12;
  assign new_n14 = ~new_n9 & ~new_n13;
  assign new_n15 = ~d & ~new_n14;
  assign new_n16 = ~a & b;
  assign new_n17 = a & ~b;
  assign new_n18 = ~new_n16 & ~new_n17;
  assign new_n19 = d & ~new_n18;
  assign new_n20 = ~c & new_n19;
  assign F0 = new_n15 | new_n20;
endmodule


