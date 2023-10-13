class NewException extends Exception {
    NewException() {
        super();
    }
    NewException(String msg, Throwable cause, boolean b1, boolean b2) {
        super(msg, cause, b1, b2);
    }
    NewException(String msg, Throwable cause) {
        super(msg, cause);
    }
    NewException(String msg) {
        super(msg);
    }
    NewException(Throwable cause) {
        super(cause);
    }
}
