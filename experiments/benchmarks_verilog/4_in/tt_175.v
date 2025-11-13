// Benchmark "ex" written by ABC on Fri Jul  4 15:14:14 2025

module top (
    a, b, c, d,
    F0  );
  input  a, b, c, d;
  output F0;
  wire new_n6, new_n7, new_n8, new_n9, new_n10, new_n11;
  assign new_n6 = ~b & c;
  assign new_n7 = b & ~c;
  assign new_n8 = ~new_n6 & ~new_n7;
  assign new_n9 = ~b & ~c;
  assign new_n10 = ~a & new_n9;
  assign new_n11 = new_n8 & ~new_n10;
  assign F0 = ~d & ~new_n11;
endmodule


