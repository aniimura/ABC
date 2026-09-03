# [Falt1] トラックのゴール(2026-09-04 設定・確定、2026-09-04 全項目確定)

対象: G. Faltings, *p-Adic Hodge Theory*, JAMS Vol.1 No.1 (1988) pp.255-299
(物理46ページ、`papers.json`に短縮タグ`Falt1`で登記済み、`0_Source`にPDF/txtあり)。

## 0. ゴール(現在地)

> **Falt1 Chapter I §1 0/2, §2 0/4, §3 0/2, §4 0/5,
> Chapter II §1 0/2, §2 0/4, §3 0/2,
> Chapter III §1 0/3, §2 0/2, §3 0/3, §4 0/1 —— 合計 0/30**
>
> ★先行あり: Chapter I §4.4(Theorem 4.4)は[LocProP]形式化の一環で
> `Interface/Falt1/PAdicHodgeTheory.lean`・`Skeleton/Falt1/ChapterI.lean`
> として既にposit済みだが、**(i)(iv)のみの縮小形**([LocProP] Lemma 2.1
> が引く部分だけ)。このトラックでの「完了」は原文の(i)-(vii)全体を
> 述べることなので、上のカウントには**含めていない**(0/5のまま)。
> §4.4に着手する際は、この既存ファイルを土台に(ii)(iii)(v)(vi)(vii)まで
> 拡張する。

## 1. 構造(2026-09-04 確定——旧版の3つの未確認点をすべて解消した)

★★**Chapter見出しは"CHAPTER"の語を伴わない裸のローマ数字**("I""II""III"
が単独行に出るだけ)。原文 物理p.4(印字p.256)の(e)段落「... In §I we study
the case of good reduction ... In §II we construct the theory 𝒳*(X) ...
Finally in §III we treat the case of bad reduction ...」で著者自身が
3章構成を明言しており、実測(下表)と完全に一致する。

★★**項目の番号付けは "N.M. Kind" の順**(番号が先、種別が後)。LocProP等の
"Kind N.M"(種別が先)と**逆**。`tools/paper-items.mjs`はそのまま使えない
(専用の抽出スクリプトが必要——再現用ロジックは本ファイル末尾のメモ参照)。

★★**pageOffsetを訂正した(254→253)**: 物理p.1はJSTORの表紙(著者・出典
スタンプのみ、本文ページ番号なし)で、物理p.2が印字p.255(論文本文の
実際の1ページ目)。150dpiで物理p.1-3を目視し、物理p.3の見出し
"256 GERD FALTINGS"で確定した。★[LocProP]形式化時に登記した
`theorem_I_4_4.src`の`pdfPage := 17`は260dpi目視による直接確認だった
ため訂正の影響を受けない(訂正後の式 印字270-253=17 で無矛盾も確認済み)。

Chapter III §4 は末尾を確認した結果、**4.1 Theoremの1件のみ**(物理p.45、
印字p.298)で、直後に p.299 の BIBLIOGRAPHY が続き論文が終わる。
→ 総項目数は **30 で確定**(旧版の「30+」の「+」は不要)。

## 2. 全項目表(物理ページは pageOffset=253 で算出・機械抽出。個別の
   260dpi目視は §4.4 以外まだ)

### Chapter I(good reduction、4節・13項目)

| 節 | 項目 | 種別 | 物理p. | 内容の頭出し(OCR、地の文の空白喪失あり) |
|---|---|---|---|---|
| §1 Ramification theory for discrete valuation rings | 1.1 | Lemma | 4 | For any extension V⊂W... the natural map... |
| | 1.2 | Theorem | 5 | If V⊂W is any extension and 𝒥'V denotes the normalization... |
| §2 Almost unramified extensions | 2.1 | Definition | 6 | Suppose A is a ring, B an A-algebra. B is called an almost... |
| | 2.2 | Theorem | 7 | Suppose B=A+mB is an almost etale covering of A, C an... |
| | 2.3 | Theorem | 7 | Suppose I⊂A is a nilpotent ideal, B an almost etale covering... |
| | 2.4 | Theorem | 8 | Suppose B is an almost etale covering of A. |
| §3 Good reduction | 3.1 | Theorem | 10 | Suppose S is a normal finite R-algebra such that S[1/p]... |
| | 3.2 | Theorem | 12 | Suppose S is a finite normal torsionfree R-algebra such that... |
| §4 Differentials and cohomology | 4.1 | Theorem | 13 | (i) The map Ω_{R/V}⊗R̄→... induces almost isomorphisms |
| | 4.2 | Theorem | 15 | (i) H^i(Δ,R̂)≅Λ^i((R⊗_V V̄)^(-1)d)⊕(rest)... |
| | 4.3 | Theorem | 17 | There exists a Γ-equivariant functorial extension |
| | **4.4** | **Theorem** | **17** | **(i) H^i(Δ,R̂)≅Λ^i((R⊗_V V̄)^(-1)d)⊕(rest)...**★先行あり(縮小形) |
| | 4.5 | Theorem | 19 | (i) The spectral sequence |

### Chapter II(𝒳*(X) の構成、3節・8項目)

| 節 | 項目 | 種別 | 物理p. | 内容の頭出し |
|---|---|---|---|---|
| §1 Construction of 𝒢* | 1.1 | Theorem | 23 | (i) There exists a spectral sequence |
| | 1.2 | Theorem | 25 | (i) There exist spectral sequences |
| §2 The isomorphism with étale cohomology | 2.1 | Lemma | 27 | Any x∈X is contained in an open U⊂X which is a K(π,1). |
| | 2.2 | Corollary | 27 | For X as above, x∈X, Spec(𝒪_{X,x̄}⊗V K̄) is a K(π,1). |
| | 2.3 | Lemma | 28 | Suppose X is smooth over V, D⊂X a divisor with normal... |
| | 2.4 | Theorem | 30 | The transformation (X proper and smooth, D⊂X a divisor... |
| §3 Relations to Hodge cohomology | 3.1 | Theorem | 35 | Suppose X is proper and smooth over V, D⊂X a divisor with... |
| | 3.2 | Theorem | 36 | Suppose X is proper and smooth over V=W(k)(Witt vectors)... |

### Chapter III(bad reduction、4節・9項目)

| 節 | 項目 | 種別 | 物理p. | 内容の頭出し |
|---|---|---|---|---|
| §1 Commutative algebra | 1.1 | Lemma | 37 | There exists an e independent of n such that the normalization... |
| | 1.2 | Proposition | 38 | Suppose R has one system of good units for some (infinite)... |
| | 1.3 | Theorem | 39 | Suppose that R is a V-algebra of finite type, geometrically ir-... |
| §2 Global methods | 2.1 | Definition | 40 | A stable punctured curve of type (g,r) over an algebraically... |
| | 2.2 | Lemma | 41 | If V is a complete discrete valuation ring with fraction field K,... |
| §3 Rigid coverings | 3.1 | Definition | 41 | A V-map f:X→Y is called a rigid étale covering if |
| | 3.2 | Theorem | 42 | Suppose f:U.→Y is a rigid étale hypercovering, a coherent... |
| | 3.3 | Theorem | 43 | Suppose R is a V-algebra of finite type, smooth in characteristic... |
| §4 Intermediate cohomology | 4.1 | Theorem | 45 | Suppose X is a smooth proper K-scheme, D⊂X a divisor with... |

## 3. [LocProP]が直接引用する箇所(確認済み・実装済み)

| LocProP側 | Falt1側 | 状態 |
|---|---|---|
| [LocProP] Lemma 2.1 | Chapter I, Theorem 4.4, (i)(iv) | ★★物理p.17(印字p.270)で260dpi目視確認済み。実装済み(縮小形) |
| [LocProP] Lemma 2.3 | Chapter I(Leray-Serre + 前2つ、Theorem 4.3?) | 未確認・未着手 |
| [LocProP] 全体(§2の土台) | almost étale extension理論(§I-2) | 未確認・未着手 |

## 4. 着手順序(葉から)

leaf-first の原則により、他の項目に依存しない・短い証明から着手する。
候補(未検証・着手時に依存グラフで裏取りすること):
1. **Chapter I §1**(1.1 Lemma, 1.2 Theorem)—— almost étale extension理論
   の手前の純代数的補題で、他の Falt1 内項目への依存が最も少なそう。
2. **Chapter III §2.1**(Definition, stable punctured curve)—— 定義のみ
   なら証明義務が無く着手しやすい。
3. Chapter I §4.4 の縮小形を(i)-(vii)全体に拡張(LocProPの土台なので
   優先度は高いが、(ii)(iii)(v)(vi)(vii)は almost étale extension の
   深い性質を要する可能性があり、葉ではない懸念がある——着手前に
   `node tools/graph.mjs`で確認)。

## 5. 番号付け規則の再現手順(scratchpad/falt1-items.mjsの後継)

`scratchpad/falt1-full-index.mjs`(本セッションで作成、再現可能):
- 章境界: 行が正規表現 `^\s*(I|II|III)\s*$` に一致する単独行を探す
  (`CHAPTER`という語は原文に無い)。
- 節見出し: `^\s*(\d+)\.\s*[A-Z][a-z]` (空白喪失により単語がくっつくが
  見出し行自体は検出できる)。
- 項目: `^\s*(\d+)\.(\d+)\.\s*(Theorem|Proposition|Lemma|Corollary|
  Definition|Remark)`("N.M. Kind"の順、LocProPと逆)。
- 印字ページ: 走り込み見出し行 `NNN GERD FALTINGS`(偶数ページ)・
  `p-ADIC HODGE THEORY NNN`(奇数ページ)を全行から収集し、各項目の
  直前・直後で挟んで判定。
- 物理ページ = 印字ページ - 253(§1で確定した値)。

## 6. 未検証事項(次の個別着手時につぶすこと)

- 上表の物理ページ番号は§4.4以外、機械抽出のみ(260dpi目視未実施)。
  逐語照合(`原文 (Falt1 p.N):`)を使う前には必ず260dpiで該当ページを
  確認すること(JSTORスキャンは地の文の単語間スペースを失うが、定理の
  数式・番号は無事——Theorem 4.4の実測で確認済みの傾向)。
- 各項目が Falt1 内の他項目にどう依存するか(依存グラフの辺)は未調査。
  着手順序の葉判定は個別に`node tools/graph.mjs --owner Falt1`で
  裏取りすること。

関連: [[measure-mathlib-before-skeleton]] / `ResearchPaper/mathlib-gap.json`の
`locprop-perfectoid-hodge-tate`(この論文が主典拠)
