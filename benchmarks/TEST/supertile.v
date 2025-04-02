module top (A, B, AND_out, OR_out);

input A;     // Input for FO1 and FO2
input B;     // Input for FO1 and FO2
output AND_out;  // Output from AND gate
output OR_out;    // Output from OR gate

// FO1 Logic
wire FO1_inv;  // Inverted signal from FO1
wire FO1_cx;   // Direct signal from FO1

// FO2 Logic
wire FO2_buf;  // Buffered signal from FO2
wire FO2_cx;   // Direct signal from FO2

// FO1 logic
assign FO1_inv = A;   // Invert input A for FO1 (INV)
assign FO1_cx = A;     // Direct connection from FO1

// FO2 logic
assign FO2_buf = B;    // Buffered input B for FO2 (BUF)
assign FO2_cx = B;     // Direct connection from FO2

// AND gate logic
assign AND_out = FO2_buf & FO1_cx;  // FO1 → CX → AND, FO2 → BUF → AND

// OR gate logic
assign OR_out = FO2_cx | FO1_inv;   // FO1 → INV → OR, FO2 → CX → OR

endmodule