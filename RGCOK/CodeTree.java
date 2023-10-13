package huffman;
import java.util.ArrayList;
import java.util.List;
public final class CodeTree {
	public final InternalNode root;  // Not null
	// Stores the code for each symbol, or null if the symbol has no code.
	// For example, if symbol 5 has code 10011, then codes.get(5) is the list [1, 0, 0, 1, 1].
	private List<List<Integer>> codes;
	// Every symbol in the tree 'root' must be strictly less than 'symbolLimit'.
	public CodeTree(InternalNode root, int symbolLimit) {
		if (root == null)
			throw new NullPointerException("Argument is null");
		this.root = root;
		codes = new ArrayList<List<Integer>>();  // Initially all null
		for (int i = 0; i < symbolLimit; i++)
			codes.add(null);
		buildCodeList(root, new ArrayList<Integer>());  // Fills 'codes' with appropriate data
	}
	private void buildCodeList(Node node, List<Integer> prefix) {
		if (node instanceof InternalNode) {
			InternalNode internalNode = (InternalNode)node;
			prefix.add(0);
			buildCodeList(internalNode.leftChild , prefix);
			prefix.remove(prefix.size() - 1);
			prefix.add(1);
			buildCodeList(internalNode.rightChild, prefix);
			prefix.remove(prefix.size() - 1);
		} else if (node instanceof Leaf) {
			Leaf leaf = (Leaf)node;
			if (leaf.symbol >= codes.size())
				throw new IllegalArgumentException("Symbol exceeds symbol limit");
			if (codes.get(leaf.symbol) != null)
				throw new IllegalArgumentException("Symbol has more than one code");
			codes.set(leaf.symbol, new ArrayList<Integer>(prefix));
			
		} else {
			throw new AssertionError("Illegal node type");
		}
	}
	public List<Integer> getCode(int symbol) {
		if (symbol < 0)
			throw new IllegalArgumentException("Illegal symbol");
		else if (codes.get(symbol) == null)
			throw new IllegalArgumentException("No code for given symbol");
		else
			return codes.get(symbol);
	}
	// Returns a string showing all the codes in this tree. The format is subject to change. Useful for debugging.
	public String toString() {
		StringBuilder sb = new StringBuilder();
		toString("", root, sb);
		return sb.toString();
	}
	private static void toString(String prefix, Node node, StringBuilder sb) {
		if (node instanceof InternalNode) {
			InternalNode internalNode = (InternalNode)node;
			toString(prefix + "0", internalNode.leftChild , sb);
			toString(prefix + "1", internalNode.rightChild, sb);
		} else if (node instanceof Leaf) {
			sb.append(String.format("Code %s: Symbol %d%n", prefix, ((Leaf)node).symbol));
		} else {
			throw new AssertionError("Illegal node type");
		}
	}
}
