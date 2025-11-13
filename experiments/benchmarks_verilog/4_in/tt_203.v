// Benchmark "ex" written by ABC on Fri Jul  4 15:14:15 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6, new_n7, new_n8, new_n9;
  assign new_n6 = ~c & d;
  assign new_n7 = ~b & new_n6;
  assign new_n8 = c & ~d;
  assign new_n9 = b & new_n8;
  assign F0 = new_n7 | new_n9;
endmodule


