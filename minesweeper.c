#include "minesweeper.h"
#include <stdio.h>
#include <stdlib.h>

void select_recursion(struct cell* game, struct cell* cellptr, int dim, int* dim_sizes);

int checkWinCondition(struct cell * game, int dim, int * dim_sizes);

void initialise_game_array(struct cell* game, int dim, int* dim_sizes);

struct cell create_cell(int i, int m, struct cell* currentCell, int totalCells);

void initialise_mines(struct cell* game, int dim, int* dim_sizes, int num_mines, int** mined_cells);

int get_num_cells(int dim, int* dim_sizes);

void initialise_adjacents(struct cell* game, int dim, int* dim_sizes);

void initialise_hints(struct cell* game, int dim, int* dim_sizes, int num_mines, int** mined_cells);

int get_cellIndex(int dim,int* dim_sizes, int* coords);


void init_game(struct cell * game, int dim, int * dim_sizes, int num_mines, int ** mined_cells) {

    // MARK: Initialise game array
    initialise_game_array(game, dim, dim_sizes);

    // MARK: Initialise mines
    initialise_mines(game, dim, dim_sizes, num_mines, mined_cells);

    // MARK: Find adjacent cells
    initialise_adjacents(game, dim, dim_sizes);

    // MARK: Initialise hints
    initialise_hints(game, dim, dim_sizes, num_mines, mined_cells);

    return;
}

void initialise_game_array(struct cell* game, int dim, int* dim_sizes){
    // emptyCellIndex keeps track of where the next empty slot to store a cell is
    int emptyIndex = 0;

    // For each dimension, we want to multiply the last dimension by dim_size
    for(int i = 0; i < dim; i++){

        // The total number of cells we have stored so far
        int totalCells = 1;
        
        // Initialise totalCells variable to the correct number
        for(int j = i-1; j >=0; j--){
            totalCells *= dim_sizes[j];
        }

        // For every cell we've stored; we want to multiply that number by the size of the next dimension
        for(int m = 0; m < dim_sizes[i]; m++){
            for(int k = 0; k < totalCells; k++){

                struct cell* currentCell = &game[(totalCells * m) + k];


                if (m == 0 && i != 0) {
                    // If m == 0, that means we are still in the first index of this dimension
                    // In this case, we don't have to create a new struct and can add to the existing cell's coordinates
                    currentCell->coords[i] = m;
                } else {
                    
                    // If we are in the first dimension, or past the first index of this dimension
                    // We create a cell and store it at the next emptyIndex
                    game[emptyIndex] = create_cell(i, m, currentCell, totalCells);

                    // And we increment the emptyIndex
                    emptyIndex++;

                }


            }


        }
        

    }

}

struct cell create_cell(int i, int m, struct cell* currentCell, int totalCells){


    struct cell newCell;

    struct cell* cellptr = &newCell;

    // Initialise all variables
    cellptr->selected = 0;

    cellptr->mined = 0;

    cellptr->num_adjacent = 0;

    cellptr->hint = 0;
    // Current index for this dimension is given by m
    cellptr->coords[i] = m;

    // Find the coordinate of the other dimensions
    for(int n = i-1; n >= 0; n--){
        // We fetch the coordinate exactly totalCells number of cells before this cell
        // The cell found at this pointer should have the same coordinates in all previous dimensions, and one decrement than the current dimension.
        cellptr->coords[n] = (currentCell-totalCells)->coords[n];
    }

    return newCell;

}

void initialise_mines(struct cell* game, int dim, int* dim_sizes, int num_mines, int** mined_cells) {

    // For each mine
    for(int n = 0; n < num_mines; n++){

        int indexOfMine = get_cellIndex(dim, dim_sizes, mined_cells[n]);

        // Change the mined variable to 1 for this cell;
        struct cell* cellptr = &game[indexOfMine];

        cellptr->mined = 1;

    }
}

int get_num_cells(int dim, int* dim_sizes){

    int num_cells = 1;
    
    for(int q = 0; q < dim ; q++){
        num_cells *= dim_sizes[q];
    }

    return num_cells;

}

void initialise_adjacents(struct cell* game, int dim, int* dim_sizes){

    // At this point all cells should have had their coordinates initialised
    int num_cells = get_num_cells(dim, dim_sizes);

    // For every cell we have
    for(int r = 0; r < num_cells; r++){

        struct cell* cell_r = &game[r];

        // We want to compare dimensions of this cell to every other cell

        for(int s = 0; s < num_cells; s++){

            struct cell* cell_s = &game[s];

            // If this is the same cell, we skip
            if(r == s){
                continue;
            }

            int adjacentCoords = 0;

            // For each dimension/coordinate that these cells have
            for(int t = 0; t < dim; t++){

                // Check if the difference between coordinate t between cell s and r, are 1 or less
                int abs_diff = abs(cell_r->coords[t] - cell_s->coords[t]);

                if(abs_diff <= 1){

                    adjacentCoords++;

                }
                

            }

            // After we've checked all the coordinates pairs between cell s and cell r
            // If the number of adjacent coordinates equal to the number of dimensions we have
            // Then the two cells are indeed adjacent

            if(adjacentCoords == dim){

                // We want to append a pointer to this cell to it's adjacent array
                cell_r->adjacent[cell_r->num_adjacent] = cell_s;

                // Increment cell_r's num_adjacent
                cell_r->num_adjacent++;

            }



        }


    }
}

void initialise_hints(struct cell* game, int dim, int* dim_sizes, int num_mines, int** mined_cells){

    int num_cells = get_num_cells(dim, dim_sizes);
    
    // For each cell
    for(int u = 0; u < num_cells; u++){

        struct cell* cell_u = &game[u];

        // For each mine
        for(int v = 0; v < num_mines; v++){

            int indexOfMine = get_cellIndex(dim, dim_sizes, mined_cells[v]);

            // Here with the index of this mine; we can now make the comparison of whether it is adjacent

            if(u == indexOfMine) {
                // The cell we are checking is the same cell as the mine
                // Skip and continue to next mine
                continue;
            }

            struct cell* cell_v = &game[indexOfMine];


            // Begin checking for adjacency 

            int adjacentCoords = 0;

            // For each dimension/coordinate that these cells have
            for(int t = 0; t < dim; t++){

                // Check if the difference between coordinate t between cell s and r, are 1 or less


                int abs_diff = abs(cell_u->coords[t] - cell_v->coords[t]);

                if(abs_diff <= 1){

                    adjacentCoords++;

                }
                
            }

            // After we've checked all the coordinates pairs between cell s and cell r
            // If the number of adjacent coordinates equal to the number of dimensions we have
            // Then the two cells are indeed adjacent

            if(adjacentCoords == dim){

                // The cell is indeed adjacent to this mine
                // Increment the hint number for this cell
                cell_u->hint++;
                

            }




        }



    }
}

int get_cellIndex(int dim,int* dim_sizes, int* coords){

    // Translate coordinate into cellIndex

    int cellIndex = 0;
    int dim_multiplier = 1;

    // Since each cell have dim dimensions
    // For each coordinate that this cell has
    for(int p = 0; p < dim; p++){
        
        // Index of the cell is given by the cumulative of each cell's coordinate multiplied by the dimensions multiplier
        cellIndex += coords[p] * dim_multiplier;

        // Update the dimension multiplier 
        dim_multiplier *= dim_sizes[p];

    }

    return cellIndex;


}

int select_cell(struct cell * game, int dim, int * dim_sizes, int * coords) {

    // Check if the coordinates are within the dimensions given
    for(int d = 0; d < dim; d++){
        if(coords[d] < 0 || coords[d] >= dim_sizes[d]){

            // Coordinate out of bounds
            return 0;

        }

    }


    // Translate coordinate into cellIndex

    int cellIndex = get_cellIndex(dim, dim_sizes, coords);

    // Create pointer to the cell
    struct cell* cellptr = &game[cellIndex];

    if(cellptr->mined == 1){
        // This cell is mined and the game is finished
        cellptr->selected = 1;
        return 1;
    }

    // If there are no mines adjacent to this cell
    if(cellptr->hint == 0){

        select_recursion(game, cellptr, dim, dim_sizes);

    }


    // If there ARE adjacents to this cell
    // Simply select it
    cellptr->selected = 1;


    // Check win condition
    if(checkWinCondition(game, dim, dim_sizes)){

        // The player has selected all cells that are not mines.
        return 2;
    } 



    return 0;
}

int checkWinCondition(struct cell * game, int dim, int * dim_sizes){

    // Check how many cells we have
    int num_cells = get_num_cells(dim, dim_sizes);


    for(int i = 0; i < num_cells; i++){

        struct cell * cellptr = &game[i];

        if(cellptr->selected == 0 && cellptr->mined == 0){
            // A cell that is not mined has not been selected
            // Player has not won

            return 0;
        }


    }
    
    // If every cell  that is not a mine, has been selected
    // Player has won
    // Returns 1 the player has won
    return 1;


}

void select_recursion(struct cell* game, struct cell* cellptr, int dim, int* dim_sizes){

    // This cell has already been selected, ignore
    if(cellptr->selected == 1){
        // printf("Returned\n");
        return;
    }

    // Next, select this
    cellptr->selected = 1;

    // Now if this cell has no adjacent mines
    if(cellptr->hint == 0){

        // Call the select_recursion function on all adjacent cells of this cell
        for(int i = 0; i < cellptr->num_adjacent; i++){

            int cellIndex = get_cellIndex(dim, dim_sizes, cellptr->adjacent[i]->coords);

            // At this point, we should have the cellIndex for one of the adjacent cells
            struct cell* newCellPtr = &game[cellIndex];

            select_recursion(game, newCellPtr, dim, dim_sizes);

        }


    }



}

