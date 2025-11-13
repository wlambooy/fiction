// Benchmark "ex" written by ABC on Fri Jul  4 15:14:11 2025

module top (
    a, b, c,
    F0  );
  input  a, b, c;
  output F0;
  wire new_n5, new_n6, new_n7, new_n8;
  assign new_n5 = b & ~c;
  assign new_n6 = a & new_n5;
  assign new_n7 = ~b & c;
  assign new_n8 = ~a & new_n7;
  assign F0 = new_n6 | new_n8;
endmodule


