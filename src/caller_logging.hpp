#pragma once
#include <iostream>

// Consolidated logging macro with unified format
#define LOG_CONTEXT(type, msg) \
    std::cout << type << " [" << __FILE__ << ":" << __LINE__ << "] " << msg << std::endl

// Helper macros for specific log types
#define LOG_CALLER(function_name, node_id, start_pos, end_pos) \
    LOG_CONTEXT("CALL", "Node=" << node_id << ", Range=[" << start_pos << "," << end_pos << ")")

#define LOG_COORDS_CALLER(node_id, position) \
    LOG_CONTEXT("COORD", "Node=" << node_id << ", Pos=" << position)

#define MARK_PROBLEMATIC_BLOCK(node_id, block_id) \
    LOG_CONTEXT("BLOCK", "Node=" << node_id << ", Block=" << block_id)
