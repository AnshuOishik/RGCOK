package huffman;
public final class InternalNode extends Node {
	public final Node leftChild;  // Not null
	public final Node rightChild;  // Not null
	public InternalNode(Node leftChild, Node rightChild) {
		if (leftChild == null || rightChild == null)
			throw new NullPointerException("Argument is null");
		this.leftChild = leftChild;
		this.rightChild = rightChild;
	}	
}
