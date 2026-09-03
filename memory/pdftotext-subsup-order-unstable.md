---
name: pdftotext-subsup-order-unstable
description: pdftotextは同じ下付き+上付き記号でも出現ごとに文字順が変わる(逐語照合の対策)
metadata:
  type: project
---

pGC §2-§4 の逐語化(2026-09-04、第 1467 前後)で実測: `pdftotext -enc UTF-8`(layout/default
いずれも同じ)は、**同じ記号 `Γ_K^v` が同一段落内で複数回出現しても、下付き・上付きの
出力順序が出現ごとに違う**——1回目は `ΓvK`(上付き v が先)、2回目・3回目は `ΓKv`
(下付き K が先)。`Γ_K^ab` はさらに別の壊れ方で `ΓaKb`(上付き "ab" が "a" と "b" に
分裂し、下付き K を挟む)になることもある。同じページの同じ論文内でも一貫しない。

**Why**: グリフの画面上の絶対座標を読み取り順に並べる pdftotext の内部アルゴリズムが、
上付き・下付きそれぞれの水平方向のオフセットに依存して順序を決めているため。フォント・
文字の高さで微妙にずれ、同じ論理記号でも版面上の実際の座標が違えば順序が変わりうる。

**How to apply**: 添字付き記号を Lean の docstring 引用(`> ...`)や `1_Structured` の
`.verbatim` に書くときは、**思い込みで「自然な」下付き→上付き順を仮定しない**。必ず
その**具体的な出現箇所**を `pdftotext -enc UTF-8`(layout/default/raw の3モード)で
実測してから書く。1_Structured の HTML では `<sup>`/`<sub>` タグの**書く順序**で
射影後の文字列順が決まる(タグは全部剥がされるので、DOM 順=文字列順)ので、
実測した順序に合わせて `<sup>`/`<sub>` を並べ替えればよい。Lean のコメントでは
`_`/`^` は単純に削除されるだけ(`leanQuoteProjection`)なので、同じ理由で
生の文字順を実測どおりに書く。

自作の検証補助: `check.mjs` と同じ `squash`/`normalize`/`matchProjection` ロジックを
複製した小さな node スクリプトを作り、`.verbatim` や Lean 引用の候補文字列を
実際の3モードのテキストに対して events で二分探索マッチさせると、
`node tools/check.mjs --structured`/`--lean` を都度フルで回すより高速に往復できる。

関連: [[gate-shell-pdftotext-differs]]、[[txt-cache-vs-live-pdftotext-differs]]
