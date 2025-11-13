// Benchmark "ex" written by ABC on Fri Jul  4 15:14:16 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6, new_n7, new_n8, new_n9, new_n10, new_n11, new_n12, new_n13,
    new_n14, new_n15, new_n16, new_n17;
  assign new_n6 = ~a & b;
  assign new_n7 = a & ~b;
  assign new_n8 = ~new_n6 & ~new_n7;
  assign new_n9 = ~c & d;
  assign new_n10 = c & ~d;
  assign new_n11 = ~new_n9 & ~new_n10;
  assign new_n12 = ~new_n8 & ~new_n11;
  assign new_n13 = a & b;
  assign new_n14 = ~a & ~b;
  assign new_n15 = ~new_n13 & ~new_n14;
  assign new_n16 = ~d & ~new_n15;
  assign new_n17 = c & new_n16;
  assign F0 = new_n12 | new_n17;
endmodule


