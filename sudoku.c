/*
 *      sudoku.c
 *
 *      Created: 08. March. 2024
 *      Author: Frank Bjørnø
 *
 * Purpose: 
 *      Solve a 4x4 sudoku using graph theory. The focus is on solving a soduko and not on 
 *      providing a good user experience. As a consequence, the sudoku problem is hard coded
 *      into the main function. This approach avoids cluttering the code document with input
 *      and error checking functions, but leaves the program ill prepared to deal with sudokus
 *      containing improper values. However, a proper sudoku solver would more than likely be
 *      graphic and interactive rather than text based and handle input and error checking 
 *      differently.
 * 
 * License:
 * 
 *          Copyright (C) 2024 Frank Bjørnø
 *
 *          1. Permission is hereby granted, free of charge, to any person obtaining a copy 
 *          of this software and associated documentation files (the "Software"), to deal 
 *          in the Software without restriction, including without limitation the rights 
 *          to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
 *          of the Software, and to permit persons to whom the Software is furnished to do 
 *          so, subject to the following conditions:
 *        
 *          2. The above copyright notice and this permission notice shall be included in all 
 *          copies or substantial portions of the Software.
 *
 *          3. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
 *          INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
 *          PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
 *          HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
 *          CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
 *          OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Description:
 *          This program uses conecepts from graph theory to solve a 4x4 sudoku. Below follows a 
 *          brief description of the basic concepts used to find a solution. Basic familiarity with 
 *          graph theory and graph coloring is assumed. The expressions cell color and cell value 
 *          is used interchangeably.
 *
 *          The logic of the program relies heavily ony an adjacency matrix and a binary 
 *          representation of the sudoku.
 *
 *          The adjacency matrix for the 4x4 sudoku is a hard coded 16 elements array of 16 bit 
 *          unsigned integers. These are not numbers in any meaningful way, just data used to 
 *          describe relationships. A cell has as its neighbors all the cells in the same row, 
 *          column and box. In addition, all the neighbors of other cells with the same color 
 *          count as neighbors - as do all cells that has already been assigned a color. In the 
 *          adjacency matrix below, row number represent cells and column numbers represent its 
 *          neighbors.
 *
 *            : 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
 *          --:-----------------------------------------------   
 *           0: 1  1  1  1  1  1  0  0  1  0  0  0  1  0  0  0 
 *           1: 1  1  1  1  1  1  0  0  0  1  0  0  0  1  0  0
 *           2: 1  1  1  1  0  0  1  1  0  0  1  0  0  0  1  0
 *           3: 1  1  1  1  0  0  1  1  0  0  0  1  0  0  0  1
 *           4: 1  1  0  0  1  1  1  1  1  0  0  0  1  0  0  0
 *           5: 1  1  0  0  1  1  1  1  0  1  0  0  0  1  0  0
 *           6: 0  0  1  1  1  1  1  1  0  0  1  0  0  0  1  0
 *           7: 0  0  1  1  1  1  1  1  0  0  0  1  0  0  0  1
 *           8: 1  0  0  0  1  0  0  0  1  1  1  1  1  1  0  0
 *           9: 0  1  0  0  0  1  0  0  1  1  1  1  1  1  0  0
 *          10: 0  0  1  0  0  0  1  0  1  1  1  1  0  0  1  1
 *          11: 0  0  0  1  0  0  0  1  1  1  1  1  0  0  1  1
 *          12: 1  0  0  0  1  0  0  0  1  1  0  0  1  1  1  1
 *          13: 0  1  0  0  0  1  0  0  1  1  0  0  1  1  1  1
 *          14: 0  0  1  0  0  0  1  0  0  0  1  1  1  1  1  1
 *          15: 0  0  0  1  0  0  0  1  0  0  1  1  1  1  1  1
 *
 *
 *          Consider the 4x4 sudoku puzzle consisting of 16 individual cells:
 *
 *          -----------
 *          |0  0|0  1|   
 *          |2  0|0  4|
 *          -----------   
 *          |0  2|0  0|   
 *          |0  4|0  0|
 *          -----------
 *
 *          A binary representation of it, where a 1 means that there is a number (other than zero) 
 *          in the corresponding cell would look like this.
 *          0  0  0  1  1  0  0  1  0  1  0  0  0  1  0  0
 *
 *          The basic idea is to combine the row of the adjacency matrix corresponding to a particular 
 *          cell - with the rows corresponding to the cells with the same color - and the binary 
 *          representation of the sudoku puzzle. Consider the cell in the second row, first column of 
 *          the sudoku puzzle with the value 2 in in, corresponding to row 4 in the adjacency matrix. 
 *          If we combine this row with the binary representation of the sudoku - and row 9 of the 
 *          adjacency matrix, corresponding to the cell in the third row, second column of the puzzle, 
 *          containing the same value as the cell under scrutiny, we get the following:
 *
 *             1  1  0  0    1  1  1  1    1  0  0  0    1  0  0  0      (row 4 in the adjacency matrix, representing the edges from cell 4)
 *          OR 0  0  0  1    1  0  0  1    0  1  0  0    0  1  0  0      (the binary representation of the sudoku matrix)
 *          OR 0  1  0  0    0  1  0  0    1  1  1  1    1  1  0  0      (row 9 of the adjacency matrix, representing the edges from cell 9)
 *          -------------------------------------------------------
 *           = 1  1  0  1    1  1  1  1    1  1  1  1    1  1  0  0	   (a)
 *
 *          This calculation is done for each cell in the function update_edges and each cell is 
 *          represented by a correponding vertex in the vertices array. The bits in the calculated 
 *          expression represent edges so that the vertex representing cell number 4 has edges to, 
 *          and is therefore neigbors with, cells 0, 1, 3, 4 (itself), 5, 6, 7, 8, 9, 10, 11, 12, 13. 
 *
 *          There are several ways to group the bits in (a). In the arrangement presented above, the 
 *          first group of four bits represents the first row, the second group represents the second 
 *          row and so forth. A different way to arrange the bits would be to take the leftmost bit 
 *          from each group of four, to represent the first column, and the rightmost bit from each 
 *          group to represent the 4th column etc. The first (upper left) box can be represented using 
 *          the two leftmost bits of the first group and the two leftmost bits of the second group etc. 
 *          These groupings are hardcoded into the masks array.
 *
 *          Any single zero in a a group of four bits identifies a cell that can be colored using the 
 *          color of the cell we are investigating. Inspecting the grouping in our example (a), it 
 *          becomes immediately clear that cell #2 (row 1, column 3) can be colored (assigned the value 2).
 *
 *          Arranging the bits in the column pattern reveals that cell 15 (row 4, column 4), would also 
 *          have been a candidate for coloring:
 *          1  1  1  1    1  1  1  1    0  1  1  0    1  1  1  0
 *
 *          And finally, arranging the bits in the box pattern would again identify cell #2 as a condidate.:
 *          1  1  1  1    0  1  1  1    1  1  1  1    1  1  0  0
 *
 *          Even though it is possible to find more than one cell that can be colored, the algorithm 
 *          only colors one cell (the first to be found) per iteration of the while loop in the 
 *          solve_sudoku function. 
 *
 *          This cell identification is handled by the function vertex_to_cell_index. We first use the 
 *          AND operation to combine the edges (a) with a mask from the masks array, for example the 
 *          mask representing the first row: 0xf000 = 0b1111 0000 0000 0000. XOR'ing with the same mask 
 *          again will set the bits corresponding to uncolored cells in the row and reset all others.
 *
 *              1  1  0  1    1  1  1  1    1  1  1  1    1  1  0  0	  (a)
 *          AND 1  1  1  1    0  0  0  0    0  0  0  0    0  0  0  0    masks[0] = 0xf000 representing the first row
 *          --------------------------------------------------------
 *            = 1  1  0  1    0  0  0  0    0  0  0  0    0  0  0  0
 *          XOR 1  1  1  1    0  0  0  0    0  0  0  0    0  0  0  0    masks[0] = 0xf000 representing the first row
 *          --------------------------------------------------------
 *            = 0  0  1  0    0  0  0  0    0  0  0  0    0  0  0  0
 *
 *          If the number of set bits (number of 1's) in this result is exactly 1, then we have 
 *          identified a group with a cell than can be colored. The number of left shifts to turn this 
 *          number into 0x8000 = 0x1000 0000 0000 0000 is also the cell number. This process is 
 *          repeated for all masks until we find a mask that yields the single 1 or until the masks 
 *          array has been exhausted, in which case we try the next vertex in the vertex array.
 *
 *          After a cell has been identified and colored, the edges of the vertices are updated with the 
 *          new information. The vertex array is then sorted and a new vertex is picked on the contidions 
 *          that it is the vertex with the highest number of edges less than 16. A vertex with 16 edges 
 *          is connected to (is neighbor to) all other cells in the sudoku matrix and can therefore not 
 *          be used to identify new cells for coloring.
 *
 *          Beyond this, the functions in this program are designed to deal with sudoku specific issues 
 *          such as missing colors, and invalid sudokus.
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>



/*
 *     strings used to make horizontal lines when printing soduko matrix
 */
uint8_t top[16] = {0xda, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xc2, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xbf, 0};
uint8_t mid[16] = {0xc3, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xc5, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xb4, 0};
uint8_t bot[16] = {0xc0, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xc1, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xc4, 0xd9, 0};



/*
 *     adjacency matrix
 */
uint16_t adj_matrix[16] = {0xfc88, 0xfc44, 0xf322, 0xf311,
                           0xcf88, 0xcf44, 0x3f22, 0x3f11,
                           0x88fc, 0x44fc, 0x22f3, 0x11f3,
                           0x88cf, 0x44cf, 0x223f, 0x113f
					 };



/*
 *     string defining rows, columns and boxes
 */
uint16_t masks[12] =      {0xf000, 0x0f00, 0x00f0, 0x000f,       //  rows
                           0x8888, 0x4444, 0x2222, 0x1111,       //  columns
                           0xcc00, 0x3300, 0x00cc, 0x0033        //  boxes
					 };
					 
/*
 *     Lookup tables for converting from index to row, column or box. This maps sudoku index to mask index
 */
uint8_t itorow[16] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
uint8_t itocol[16] = {4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7};
uint8_t itobox[16] = {8, 8, 9, 9, 8, 8, 9, 9, 10, 10, 11, 11, 10, 10, 11, 11};



/*
 *     Misc. codes used when solving the Sudoku
 */
enum return_codes{OK = 0x00, 
                  NO_COLOR      = 0x00, 
			   TRUE          = 0x01, 
			   INVALID       = 0x01, 			   
			   FEW_CLUES     = 0x02, 			   
			   MAX_COUNT     = 0x10,
			   NOT_EXCLUSIVE = 0x10
			  };


/*
 *     An array of vertex structures makes it easier to control data.
 */
typedef struct
{
	uint16_t _cell, _color, _edges, _count;	
} Vertex;

enum {NVERTICES = 0x10};
Vertex vertices[NVERTICES];



/*
 *     Count the number of set bits in an unsigned 16-bit integer
 *
 *     \param word     An unsigned 16-bit integer
 *
 *     Example:        word = 0xfda8 = 0b1111 1101 1010 1000, so bit count is 10 because there are 10 set bits (and 6 zeros).
 */
uint8_t count_bits(uint16_t word)
{
	uint8_t count = 0;
	for (uint16_t mask = 1; mask; mask <<= 1) count += ((word & mask) != 0);
	return count;
}



/*
 *     Prints a 4x4 sudoku
 *
 *     \param *s   Pointer to a 16-element array of unsigned 8-bit integers or unsigned char (uint8_t), 
 *                 representing a 4x4 sudoku.
 */
void print_sudoku(uint8_t *s)
{			
	printf("    %s\n", top);	
	printf("    %c %u  %u %c %u  %u %c\n", 0xb3,  s[0],  s[1],  0xb3,  s[2],  s[3], 0xb3);
	printf("    %c %u  %u %c %u  %u %c\n", 0xb3,  s[4],  s[5],  0xb3,  s[6],  s[7], 0xb3);
	printf("    %s\n", mid);
	printf("    %c %u  %u %c %u  %u %c\n", 0xb3,  s[8],  s[9],  0xb3, s[10], s[11], 0xb3);
	printf("    %c %u  %u %c %u  %u %c\n", 0xb3, s[12], s[13],  0xb3, s[14], s[15], 0xb3);
	printf("    %s\n", bot);	
}



/*
 *     Checks rows, columns and boxes to ensure that there are exactly 4 different colors in each row, column and box.
 *
 *     \param *s   Pointer to a 16-element array of unsigned 8-bit integers or unsigned char (uint8_t), 
 *                 representing a 4x4 sudoku.
 * 
 *     \return     An 8-bit error code where the 4 most significant bits are a number from 1 to 4, and the 4 least significant
 *                 bits are a number from 0 to 2 indicating row, column or box, respectively. So the error code 0x31 indicates
 *                 an error in column 3.
 *
 *     Note:       This function is designed to validate a finsished Sudoku only, so passing an unsolved Sudoku to this function
 *                 will give no insight whatsover into the solvability of a Sudoku.
 */
uint8_t validate_sudoku(uint8_t *s)
{	
		//  check 4 rows, 4 columns and 4 boxes 
	for (uint8_t i = 0; i < 12; i++)
	{
		uint16_t mask = masks[i];
		uint8_t table = 0;			//  tracks colors
		
			//  index refers to the Sudoku array
		for (int8_t index = 0; mask; index++)
		{
			if (mask & 0x8000 && s[index] != 0) table |= (1 << s[index]);
			mask <<= 1;	
		}
		if (count_bits(table) != 4) return (i % 4 + 1) << 4 | (i / 4);		
	}
	return OK;
}




/*
 *     Populate the vertice array
 *
 *     \param *s   Pointer to a 16-element array of unsigned 8-bit integers or unsigned char (uint8_t), 
 *                 representing a 4x4 sudoku.
 */ 
void prepare_vertices(uint8_t *s)
{
	for (uint8_t i = 0; i < 16; i++)
	{
		vertices[i]._cell = i;
		vertices[i]._color = s[i];
		vertices[i]._edges = adj_matrix[i];
	}
}



/*
 *     Compare two vertices by the _count member. Only called from qsort. By returning
 *     b._count - a._count the elements will be sorted in descending order.
 *
 *     \param *a  Pointer to the first element to be compared.
 *     \param *b  Pointer to the second element to be compared.
 *
 *     \return    b._count - a._count to sort in descending order
 */
int cmp_vertices(const void *a, const void *b)
{
	return (*(Vertex*)b)._count - (*(Vertex*)a)._count;	
}



/*
 *     Sort vertices in descending order using the cmp_vertices function defined above
 */
void sort_vertices()
{
	qsort(vertices, 16, sizeof(Vertex), cmp_vertices);
}



/*
 *     Update connections between the vertices and its neighbors
 *
 *     \param *s   Pointer to a 16-element array of unsigned 8-bit integers or unsigned char (uint8_t), 
 *                 representing a 4x4 sudoku.
 */
void update_edges(uint8_t *s)
{
		//  Make a binary representation of the matrix. If there is a number in a given cell of the sudoku, 
		//  the corresponding bit in the cells variable is set
	uint16_t cells = 0, mask = 0x8000;
	for (uint8_t *ptr = s; ptr < s + 16; mask >>= 1) cells |= mask * (*ptr++ != 0);	
	
		//  Connect each vertex to its neighbors
	Vertex *v = vertices;
	for (int c = 0; c < 16; c++, v++)
	{
			//  skip vertices without color
		if (v -> _color == 0) continue;
		
			//  add already colored vertices to _edges
		v -> _edges |= cells;
		
			//  share edges with vertices of the same color
		for (int i = 0; i < 16; i++) v -> _edges |= vertices[i]._edges * (vertices[i]._color == v -> _color);
		
			//  update _count to reflect changes in edges
		v -> _count = count_bits(v -> _edges);
	}
}



/*
 *     Pick the vertex with the highest _count that is less than 16, starting from a given vertex. The
 *     array of vertices must be sorted in descending order.
 *
 *     \param *start   The starting vertex.
 *
 *     \return         The vertex with the highest count less than 16 or NULL if either all vertices
 *                     have _count = 16 or the array of vertices has been exhausted.
 */
Vertex* pick_vertex(Vertex *start)
{		
	if (start == NULL)
	{
		start = vertices;
		for (uint8_t i = 0; i < NVERTICES && start -> _count == MAX_COUNT; i++) start++;
	}
	else
	{
		start++;
	}
	
	if (start - vertices == NVERTICES) return NULL;
	
	return start;
}



/*
 *     Look up a vertex corresponding to a cell.
 *
 *     \param cell   A cell number (0 - 15)
 *
 *     \return       The corresponding vertex, that is, the vertex whose _cell member matches the cell parameter
 *
 *     Warning:   The function does no error checking and is only called from the solve_sudoku function when
 *                the cell argument is guaranteed to be limited to the interval [0, 15].
 */
Vertex* cell_index_to_vertex(uint8_t cell)
{
	Vertex *v = vertices;
	for (uint8_t i = 0; i < NVERTICES; i++)
	{
		if (v -> _cell == cell) break;
		v++;
	}
	return v;
}


/*
 *     Calculates the position of an exclusively empty cell based on vertex edges
 *
 *     \param edges   v._edges of the vertex in question
 *
 *     \return        The position (0 - 15) of the empty cell, 16 if not exclusive or no empty cell
 */
uint8_t vertex_to_cell_index(Vertex *v)
{
	uint16_t edges = v -> _edges;
	if (edges == MAX_COUNT) return NOT_EXCLUSIVE;
	
	uint16_t position = 0;
	
		//  set all bits that represent empty cells
	for (uint8_t i = 0; i < 12; i++)
	{
		position = (edges & masks[i]) ^ masks[i];
		if (count_bits(position) == 1) break;
	}
	
		//  if more than one bit is set, the cell is not the only empty sell in the set
	if (count_bits(position) != 1) return NOT_EXCLUSIVE;
	
		//  calculate cell position
	uint8_t cell = 0;
	for (; position != 0x8000; position <<= 1) cell++;
	
	return cell;
}



/*
 *     Locate the position of a cell based on a given mask
 *
 *     \param  mask   One of the masks from the array masks
 *     \param *s      Pointer to a 16-element array of unsigned 8-bit integers or unsigned char (uint8_t), 
 *                    representing a 4x4 sudoku.  
 *    
 *     \return        The position (0 - 15) of the empty cell.
 *
 *     Note:   This function is only called from the function loose_cell, at a point where it has been 
 *             established that there is a loose cell.
 */
uint8_t mask_to_cell_index(uint16_t mask, uint8_t *s)
{	
	for (uint8_t cell = 0; cell < 16; cell++)
	{
		if (mask & 0x8000 && s[cell] == 0) return cell;
		mask <<= 1;
	}
	return NOT_EXCLUSIVE;
}



/*
 *     Calculates the mask index of a "loose cell", that is a cell, whose color can only be deduced
 *     by the fact that all other cells in the set (row, column or box) are filled. F.ex. if a row
 *     have the clues 0 3 1 2, we can deduce that the missing clue, 0, must be 4. The existence of a
 *     loose cell is uncovered by trying all masks to se if there is a row, column or box with exactly
 *     three colors.
 *
 *     \param *s   Pointer to a 16-element array of unsigned 8-bit integers or unsigned char (uint8_t), 
 *                 representing a 4x4 sudoku.
 *
 *     \return     The mask corresponding to the loose cell
 */ 
uint16_t has_loose_cell(uint8_t *s)
{
	for (uint8_t i = 0; i < 12; i++)
	{
		uint16_t mask = masks[i];
		uint8_t  colors = 0;
				
		for (uint8_t index = 0; mask; index++)
		{				
			colors |= (mask & 0x8000 && s[index] != 0) << s[index];
			mask <<= 1;
		}		
		if (count_bits(colors) == 3) return masks[i];				//  i is the index of the mask used to uncover the loose cell
	}
	return NOT_EXCLUSIVE;
}



/*
 *     Query the missing color in a set (row, column or box) defined by cell index.
 *     If the sudoku is valid, any missing color will be the same in the cell's row, 
 *     column and box. Otherwise, the sudoku is invalid
 *
 *     \param  ci  Cell index, that is index of a loose cell.
 *     \param *s   Pointer to a 16-element array of unsigned 8-bit integers or unsigned char (uint8_t), 
 *                 representing a 4x4 sudoku.
 *
 *     \return     The missing color (1 - 4) or NO_COLOR = 0 if the color is not consistently missing
 *                 from row, column and box.
 *
 *     Note:   There is a hypothetical possibility that more than one color could be missing. If so, the function 
 *             would return the lowest number. The issue could be addressed by counting the bits in result before 
 *             the initialization. However, the function is only called from a point in the solve_sudoku algorithm
 *             where the existence of a loose cell has been established, meaning that there are exactly one color
 *             missing or none (in which case the sudoku is invalid).
 */
uint8_t missing_color(uint8_t ci, uint8_t *s)
{
		//  Row Mask, Column Mask and Box Mask
	uint16_t rcb[3] = {masks[itorow[ci]], masks[itocol[ci]], masks[itobox[ci]]};
	
		//  Color tables for rowm column and box
	uint8_t rc = 0, cc = 0, bc = 0;
	for (uint16_t mask = 0x8000, i = 0; mask; mask >>= 1, i++)
	{
		rc |= ((rcb[0] & mask) != 0) << (s[i] - 1);
		cc |= ((rcb[1] & mask) != 0) << (s[i] - 1);
		bc |= ((rcb[2] & mask) != 0) << (s[i] - 1);		
	}

		//  and together the inverses of the color tables
	uint8_t result = ~rc & ~cc & ~bc;
		
		//  isolate low nybble and count shifts until bit1 is set. This number +1 represents color
	if (result & 0x0f)
	{
		for (uint8_t c = 0; c < 4; c++, result >>= 1)
		{
			if (result & 1) return ++c;
		}
	}
	
	return NO_COLOR;
}



/*
 *     Solves the sudoku
 *
 *     \param *s   Pointer to a 16-element array of unsigned 8-bit integers or unsigned char (uint8_t), 
 *                 representing a 4x4 sudoku.
 */
uint8_t solve_sudoku(uint8_t *s)
{
	Vertex *v = NULL;
	uint8_t cell, color;
	
	while (TRUE)
	{
		update_edges(s);
		sort_vertices();
		
PICK_VERTEX:
		if ((v = pick_vertex(v)) == NULL) return validate_sudoku(s);	//  either OK or INVALID
		
		if (v -> _count != 0)
		{
			cell = vertex_to_cell_index(v);
			if (cell == NOT_EXCLUSIVE) goto PICK_VERTEX;			//  try next vertex
			
			color = v -> _color;			
		}
		else
		{			
			uint16_t mask = has_loose_cell(s);
			if (mask == NOT_EXCLUSIVE) return FEW_CLUES;
			
			cell = mask_to_cell_index(mask, s);
			color = missing_color(cell, s);
			if (color == NO_COLOR) return INVALID;
		}
		
		s[cell] = color;
		cell_index_to_vertex(cell) -> _color = color;
		v = NULL;
	}
	
	return OK;									//  unreachable
}



int main(void)
{	
		//  starting matrix
	uint8_t sudoku[16] = {0, 1, 0, 4,   0, 0, 1, 0,   0, 4, 0, 3,   0, 0, 0, 0};
	
		//  prepare vertices	
	prepare_vertices(sudoku);
	
		//  present initial problem
	print_sudoku(sudoku);
	
		//  attempt to solve sudoku, report any error, and present final result
	int result = solve_sudoku(sudoku);
	
	switch(result)
	{
		case OK:
			break;		
		case FEW_CLUES:
			puts("Not enough clues to solve the sudoku.");
			break;
		default:
			puts("Invalid sudoku.");
			break;
	}
	
	print_sudoku(sudoku);	
	
	
	return 0;
}
