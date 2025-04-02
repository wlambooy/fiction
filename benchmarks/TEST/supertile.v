module top (A, B, AND_out, OR_out);

input A;     // Input for FO1 and FO2
input B;     // Input for FO1 and FO2
output AND_out;  // Output from AND gate
output OR_out;    // Output from OR gate

// FO1 Logic
wire FO1_buf;  // Buffered signal from FO1
wire FO1_cx;   // Buffered signal from FO1

// FO2 Logic
wire FO2_cx;   // Buffered signal from FO2
wire FO2_inv;  // Inverted signal from FO2

// FO1 logic
assign FO1_buf = A;   // Buffer input A for FO1 (BUF)
assign FO1_cx = A;     // Buffered connection from FO1

// FO2 logic
assign FO2_cx = B;     // Buffered connection from FO2
assign FO2_inv = B;    // Invert input B for FO2 (INV)

// AND gate logic
assign AND_out = FO1_buf & FO2_cx;  // FO1 → CX → AND, FO2 → BUF → AND

// OR gate logic
assign OR_out = FO1_cx | FO2_inv;   // FO2 → INV → OR, FO2 → CX → OR

endmodule