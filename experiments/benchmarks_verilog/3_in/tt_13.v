// Benchmark "ex" written by ABC on Fri Jul  4 15:14:11 2025

module top (
    a, b, c,
    F0  );
  input  a, b, c;
  output F0;
  wire new_n5, new_n6, new_n7, new_n8, new_n9;
  assign new_n5 = ~a & b;
  assign new_n6 = ~a & ~new_n5;
  assign new_n7 = ~c & ~new_n6;
  assign new_n8 = ~b & c;
  assign new_n9 = ~a & new_n8;
  assign F0 = new_n7 | new_n9;
endmodule


