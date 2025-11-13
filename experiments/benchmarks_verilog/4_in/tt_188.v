// Benchmark "ex" written by ABC on Fri Jul  4 15:14:14 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6, new_n7, new_n8;
  assign new_n6 = ~b & ~c;
  assign new_n7 = c & ~d;
  assign new_n8 = b & new_n7;
  assign F0 = new_n6 | new_n8;
endmodule


