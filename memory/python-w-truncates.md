---
name: python-w-truncates
description: Python の io.open(path,'w') は書き込み前にファイルを 0 バイトへ切り詰める——encode で落ちると原本が消える
metadata:
  type: feedback
---

`io.open(p,'w',encoding='utf-8').write(t)` は **open した時点でファイルを 0 バイトにする**。
その後 `write` が `UnicodeEncodeError`（サロゲート対など）で落ちると、**元の中身は失われる**。

**Why:** 2026-08-21、`ABC3/Probe.lean` をこれで消した。しかも
**空ファイルは `lake env lean` が無言で通す**ので、その後の 4 回のビルドが
「成功」に見えていた。★消えたことにも、検証が空回りしていたことにも気づかなかった。

**How to apply:**
* 一時ファイルへ書いて `os.replace` する（`ResearchPaper/frdi-decomposition.json` の更新でやっている形）。
  `io.open(p+'.tmp','wb').write(s.encode('utf-8')); os.replace(p+'.tmp', p)`
* 星型外の文字（𝟙 𝒞 𝔽 𝓞 …）は Python 文字列リテラルに `\uXXXX` で書かない。`chr(0x1D7D9)` を使う。
  → [[surrogate-pair-encoding]]
* ★**Lean の検証は「出力が無いこと」ではなく「宣言が実在すること」で確かめる**。
  空ファイルもエラーを出さない。`wc -l` か `grep -c '^theorem\|^def'` を併せて見ること。
