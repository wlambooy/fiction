// Benchmark "ex" written by ABC on Fri Jul  4 15:14:13 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6, new_n7, new_n8, new_n9, new_n10, new_n11, new_n12;
  assign new_n6 = c & ~d;
  assign new_n7 = c & ~new_n6;
  assign new_n8 = ~b & ~new_n7;
  assign new_n9 = ~a & new_n8;
  assign new_n10 = a & b;
  assign new_n11 = ~c & ~d;
  assign new_n12 = new_n10 & new_n11;
  assign F0 = new_n9 | new_n12;
endmodule


