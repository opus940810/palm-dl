"""Extract candidate site mentions from whatever full texts have been cached."""
from pathlib import Path

from lit_mining import OUT, extract_texts


def main() -> None:
    files = sorted((OUT / "xml").glob("*.xml"))
    if not files:
        raise SystemExit("no cached full texts")
    texts = [p.read_text(encoding="utf-8", errors="ignore") for p in files]
    df = extract_texts(texts, [p.stem for p in files])
    df.to_csv(OUT / "candidates.csv", index=False)
    print(f"{len(files)} full texts -> {len(df)} mentions, {df.token.nunique()} tokens")


if __name__ == "__main__":
    main()
