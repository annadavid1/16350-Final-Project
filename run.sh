g++ planner.cpp -o planner.out \
-I/opt/homebrew/opt/sfml@2/include \
-L/opt/homebrew/opt/sfml@2/lib \
-lsfml-graphics -lsfml-window -lsfml-system \
-pthread

./planner.out $1