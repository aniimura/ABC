# [GenEll] トラックの長期ゴール（2026-08-16 設定）

対象: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]（物理 25 ページ、番号付き項目 **26 件**）。

本書は `PLAN.md` §0 の北極星（IUT 全体 + abc）に対する**下位の北極星**を 1 つ設定する。
`layer0-triage.md` が「★★いま着手できるもの」の 5 番に挙げた `[GenEll] Lemma 3.1`
（純古典・IUT 語彙ゼロ）が**何に向かう最下段なのか**を測って決めた。

**★すべて `.txt`（pdftotext）由来であり、目視していない。**事実2により項目名・逐語は壊れうる。
本書の数はすべて「桁を見るための数」であり、逐語には使えない。着手時に PDF 目視が要る。

---

## 0. 一行

> **★北極星（B トラック）: `[GenEll] Corollary 4.4` まで、IUTchIV が要る範囲の [GenEll] を `Found/` に `sorry` 無しで載せる。**
>
> **★同じ土台の上に立つ第二の頂: `[GenEll] Theorem 2.1`** ——「abc ⟺ ℙ¹∖{0,1,∞} 版 abc」という**古典的同値性**。**IUT を一切使わない。**

---

## 1. なぜこの 2 つが頂になるか（実測）

### 1-1. Corollary 4.4 の側 —— IUT の**入力が空虚でない**ことを言う唯一の道

[IUTchIV] Corollary 2.2 の証明（`.txt` l.3055–3100）は、与えられた楕円曲線から
初期 Θ-データを構成するために (P1)–(P7) を順に立てる。その中の (P6):

> (P6) the image of the outer homomorphism Gal(Q/F) →GL2(Fl) determined
> by the l-torsion points of EF contains the subgroup SL2(Fl) ⊆GL2(Fl).
>
> Indeed, since, by (P5), EF has bad multiplicative reduction at some valuation ∈
> Vbad_mod ≠ ∅, (P6) follows formally from (P2), (P4), and **[GenEll], Lemma 3.1, (iii)**
> [cf. the proof of the final portion of [GenEll], Theorem 3.8].

そして (P1)–(P7) が揃うことで

> there exist data “CK”, “V”, and “ϵ” such that all of the conditions of
> **[IUTchI], Definition 3.1, (a), (b), (c), (d), (e), (f), are satisfied**

★**これは `PLAN.md` §4 の G2（非空虚 witness）の穴に、この一点で正面から当たる。**
G2 は「退化 witness でも通る」ことが実証済みで、内容を保証しない。
本プロジェクトの `Skeleton/IUTchI/InitialThetaData.lean` は現在 **Definition 3.1 の条件 (a) だけ**の転写であり、
その docstring 自身が「これは initial Θ-data ではない」と書いている。
**初期 Θ-データの「本物の witness」を作る道は、(P1)–(P7) を通る以外に無く、その (P6) が Lemma 3.1 (iii) である。**

### 1-2. Theorem 2.1 の側 —— abc への最終段の**半分そのもの**で、しかも IUT 非依存

[IUTchIV] Corollary 2.3（= Diophantine Inequalities。`dependency-scale.md` が「ABC 予想を導くもの」と記録している系）の証明（`.txt` l.3608–3640）:

> One verifies immediately that the content of the statement of Corollary 2.3
> coincides precisely with the content of **[GenEll], Theorem 2.1, (i)**. Thus, it follows
> from **the equivalence of [GenEll], Theorem 2.1**, that, in order to complete the proof
> of Corollary 2.3, it suffices to verify that **[GenEll], Theorem 2.1, (ii)**, holds.

すなわち abc への最終段は 2 枚に分かれる:

| | 内容 | 誰が示すか |
|---|---|---|
| (A) | 一般の双曲曲線の Vojta 不等式 ⟺ ℙ¹∖{0,1,∞} の compactly bounded subset 版 | **[GenEll] Theorem 2.1**（Elkies / Moret-Bailly / Szpiro の技法 + noncritical Belyi maps）|
| (B) | その ℙ¹ 版を実際に示す | **IUT 本体**（Cor 3.12 → IUTchIV Cor 2.2）|

**本プロジェクトはこれまで (B) 側しか見ていない。(A) は独立に形式化でき、IUT の当否に依存しない。**
`idea2.md` ①「解けなくても数学的成果が残る問題を選ぶ」の理想形であり、
しかも (B) が完成したときに**接続がそのまま閉じる**。

---

## 2. 規模（実測。★器具を作り直して測った）

### ★2-0. 訂正の記録 —— `tools/intra-graph.mjs` の近似が結論を変えていた

`intra-graph.mjs` は項目の本体を「**次の宣言行まで**」と近似する。
節の最後の項目では**次の節の導入文を巻き込む**。実害を実測した:

| | Theorem 2.1 の到達 | 内訳 |
|---|---:|---|
| `intra-graph.mjs`（次の宣言行まで） | **17 件** | §3 の項目（Theorem 3.8 / Lemma 3.1 / …）を含む |
| 節見出しでも切る（本書） | **9 件** | **§1 と自身のみ** |

原因は §2 が Theorem 2.1 **1 項目しか持たない**こと。その本体が `Section 3:` の導入文
（`… is “rather large” [cf. Theorem 3.8]`）まで伸び、そこにある参照を Theorem 2.1 の依存として数えていた。

> ★**「Theorem 2.1 は Theorem 3.8 に依存する」は誤りである。**
> Theorem 2.1 の本体（118 行）が引く論文内項目は
> **Example 1.3 ×2 / Proposition 1.7 ×3 / Proposition 1.4 ×2 / Proposition 1.6 ×1** の 4 件だけ、
> 外部文献は **[Mzk1] ×2** だけである。
> これは規模の見積りを **17 → 9** に、しかも**最重量部（§3 の Galois 表現）を圏外に**動かす差である。

再現: `scratchpad/genell-closure.mjs`（節見出しでも本体を切る版）。
★`tools/intra-graph.mjs` 自体は**直していない**（別セッションが稼働中の共有ファイルのため）。
直すなら「本体の終端 = 次の宣言行 **または** 次の `Section N:` 見出しのうち先に来る方」。

### 2-1. 到達（節見出しで切った版）

| 根 | 到達 | 深さ | 内訳 |
|---|---:|---:|---|
| **Theorem 2.1** | **9** | 3 | §1 の 8 件 + 自身。**§3・§4 を含まない** |
| **Theorem 3.8** | 13 | 3 | §1 の 4 件 + §3 の 8 件 + 自身 |
| **Corollary 4.4** | **18** | 3 | §1 の 4 件 + §3 の 9 件 + §4 の 4 件 + 自身 |
| **Lemma 3.1** | **1** | 0 | ★**真の葉**（出次数 0） |

比較（`dependency-scale.md`）: [FrdI] は 83 項目・到達 45・辺 184・平均出次数 4.1。
**[GenEll] は項目数で FrdI の 1/3、到達で 18/26。**

### 2-2. IUTchIV が [GenEll] を直接引く回数（実測、14 種）

```
Theorem 2.1 ×10 / Example 1.3 ×7 / Proposition 1.4 ×5 / Lemma 3.5 ×5 / Definition 1.2 ×5
Proposition 3.4 ×3 / Theorem 3.8 ×2 / Definition 1.5 ×2 / Corollary 4.4 ×2
Remark 4.1.1 / Lemma 3.7 / Lemma 3.1 / Definition 3.3 / Definition 1.1 各 ×1
```

★**§1（高さ）は両方の頂の共通の土台である。** Theorem 2.1 の到達 9 件のうち 8 件が §1、
Corollary 4.4 の到達 18 件のうち 4 件が §1。

---

## 3. ★証明側の依存（statement には現れない）

`PLAN.md` §8-1 の教訓（pGC §1 の選定の誤り = *statement が言及する対象*だけを見て*証明の依存*を見なかった）を、ここで先に払う。

### 3-1. 外部文献（GenEll 本文中の出現位置つき）

| 文献 | 使う箇所 | `0_Source` |
|---|---|---|
| **[Mzk1]** Noncritical Belyi Maps, Thm 2.5 | **Theorem 2.1 の証明**（l.627） | ★**ある**（`Noncritical Belyi Maps.pdf`）|
| **[Serre]** Abelian l-adic Repr., Ch IV §3.4 Lemma 3 | **Lemma 3.1 (iv)**（l.798） | 無い（書籍 1968）|
| [FC] Degenerations of Abelian Varieties, I 2.7 / III 7.3 / V 4.5 | §3 の算術側（l.823, 928, 955）| 無い |
| [Falt] / [Merel] / [Silv1] / [Szp] / [Elkies] / [vF] | §2・§3 の背景と技法 | 無い |

★**`[EL]` は参考文献ではない。** 本文で 59 回現れ外部文献の最多に見えるが、
**楕円曲線 `E_L` の下付き文字が pdftotext で `[EL]` になったもの**である（`[EL] ∈Mell(Q)`）。
GenEll の参考文献欄に `[EL]` は無い。事実2の実例をもう 1 つ。
**外部文献を角括弧パターンで自動計数すると誤る**（`tools/bibmap.mjs` の入力にする際は要注意）。

### 3-2. mathlib（探索範囲つき。`check.mjs` A2 の S1–S4）

| 我々が要るもの | 原文側の呼称 | mathlib | 探索範囲 |
|---|---|---|---|
| 高さ | “height” | ★**ある**。`Mathlib/NumberTheory/Height/` **5 ファイル 2006 行**（`mulHeight`/`logHeight`/`Northcott`/`NumberField`/`Projectivization`）| `find Mathlib -ipath "*eight*" -name "*.lean"` |
| Sylow 論 | “l-Sylow subgroup” | ある（272 件）| mathlib 全体 grep |
| 交換子群 | “commutator subgroup” | ある（591 件）| 同上 |
| transvection | （原文は行列で書く）| ある（288 件）。`diagonal_transvection_induction_of_det_ne_zero`（体上）まで | `Mathlib/LinearAlgebra/Matrix/Transvection.lean` |
| 群の完全性 | “no nontrivial abelian quotients” | ★**無い**。`perfect` の hit は**すべて perfect matching**（グラフ理論）| mathlib 全体 case-insensitive |
| SL₂(F_l) の生成 | “generated by the matrices α, β” | ★**無い** | `LinearAlgebra/SpecialLinearGroup.lean` / `LinearAlgebra/Matrix/SpecialLinearGroup.lean` / `NumberTheory/Modular.lean` / `NumberTheory/ModularForms/*.lean` |
| 楕円曲線の l-捩れ Galois 加群 | “the Galois representation Gal(Q/L) → GL2(Zl)” | ★**無い**（`torsion` の hit 8 件はすべて 2-torsion polynomial）| `Mathlib/AlgebraicGeometry/EllipticCurve/*.lean` |
| Belyi 写像 | “Belyi map” | ★**無い**（0 件）| mathlib 全体 case-insensitive |

★`perfect` の件は「**語の一致は当てにならない**」（`layer0-triage.md` の `sharp` と同型）の実例をもう 1 つ増やす。

---

## 4. 段階 —— 各段が単独で成果になるように切る

| 段 | 内容 | 前提 | 止まってもここまでは残る |
|---|---|---|---|
| **0** | **Lemma 3.1 (i)(ii)(iii)** —— SL₂(F_l) の生成 / 完全性 / Sylow 判定 | mathlib の Sylow・commutator・transvection | ★**IUTchIV が (P6) で直接使うのは (iii) である**。mathlib へ PR できる純群論 |
| **1** | **Lemma 3.1 (iv)** —— 閉部分群 J ⊆ GL₂(ℤ_l) の持ち上げ | [Serre] Ch IV §3.4 Lemma 3 は `0_Source` に無い → **自分で証明する**（Frattini 型。mathlib の `Frattini` は 15 件）| ★これで **Theorem 3.8 の群論側が完成**する（原文の証明が使う群論は Lemma 3.1 (iv) だけ、実測）|
| **2** | **§1（高さ）** —— Definition 1.2 / Example 1.3 / Proposition 1.4 / Definition 1.5 / Proposition 1.6 / 1.7 | mathlib の Height 2006 行の上に **BD-class**（bounded discrepancy）と **compactly bounded subset** を載せる | ★**両方の頂の共通の土台**。IUTchIV が最も多く引く群（Ex 1.3 ×7 / Prop 1.4 ×5 / Def 1.2 ×5）|
| **3** | **§3 の算術側** —— Lemma 3.2 / Def 3.3 / Prop 3.4 / Lemma 3.5・3.6・3.7 | ★**l-捩れ Galois 表現が mathlib に無い**（`layer0-triage.md` の壁の種類 4 = 最も重い）。Tate 曲線・semi-stable reduction・[FC] 3 箇所 | person-month 級。★**stop-loss の主対象** |
| **4** | **Theorem 3.8 → §4 → Corollary 4.4** | 段 1 + 段 3 | 北極星（B）到達 |
| **★別線** | **Theorem 2.1** | **段 2 + [Mzk1] の転写**（★段 3 は要らない、実測）| ★**abc ⟺ ℙ¹版 abc**。IUT が偽でも真であり続ける古典的同値性 |

★**段 2 が終わった時点で、Theorem 2.1 に残るのは [Mzk1] だけになる。**
論文内の依存は §1 で尽きる（§2-0 の訂正）。[Mzk1] は `0_Source` にある。
**したがって「§3 の壁（l-捩れ Galois 表現）に当たる前に、abc 側の頂が 1 つ射程に入る。」**

---

## 5. 受理条件（既存のゲートに合わせる）

**着手前の準備（G1 の前提）**

- `ResearchPaper/papers.json` に **`GenEll` を登記**する。★現在の登記は 8 本（pGC / EtTh / IUTchIII / AbsTopIII / SemiAnbd / IUTchI / IUTchII / FrdI）で、**GenEll は未登記**。登記なしに `.src` を書くと G1 が落ちる。
- `notationRisk` は**未測定**。最初の目視（物理 p.13–14 = Lemma 3.1）で決める。★GenEll は行列を多用するので、`pdftotext` の 2 次元レイアウト破壊が EtTh/IUTchIII 型の下線問題とは**別種の危険**として出る可能性がある（実例: Lemma 3.1 の α, β, γ は `.txt` 上で行列の中身がばらける）。
- 段 4 以降に進むなら `[Mzk1]`（Noncritical Belyi Maps）も `papers.json` と `tools/_fileof.json` に足す（`_fileof.json` に現在**無い**）。

**各項目**

- `Found/` は `sorry` ゼロ・`axiom` ゼロ（G4）。止まったら**理由を定理か docstring として残す**（`sorry` を置かない）。
- `.src` は「**その原典項目を完全に実装した**」という 2 値の主張。条つき（`, (i)`）は数に入れない。
- 完成宣言の前に、**文脈を持たない子による監査**（`memory/challenger-audit-without-context.md`）。
- 進捗は機械で数える: **`tools/genell-progress.mjs`**（分母 = **IUTchIV からの需要の推移閉包**）。
  ★`tools/frdi-progress.mjs` は FrdI 決め打ちなので、共通化ではなく**別ファイル**として作る
  （FrdI 側が稼働中の共有ファイルを触らないため。共通化はどちらも止まっているときに）。

---

## 6. stop-loss

`PLAN.md` §8-3 に従い、**時間ではなく進捗の質**で切る。

- ★**段 3 が唯一「壁の種類 4（語彙の不在）」に当たる段**である。ここに入る前に段 0–2 を完成させ、
  「何が要るか」を `ProofObligation` として書き切っておく。
- 段 3 に入ってから「同じエラーに 3 回戻る／迂回路を 2 つ試して全部同じ壁」になったら、
  そこで **Theorem 2.1 の別線（段 2 + [Mzk1]）へ移る**。★この goal は頂が 2 つあるので、
  片方が詰まっても撤退にならない。

---

## 7. この目標が外れうる点（先に書く）

1. **§1 の 2 つの引用は `.txt` 由来で目視していない。** [IUTchIV] p.46 の (P6) と p.55 の Cor 2.3 の証明は、
   **目視で確かめるまで暫定**である。特に (P6) の `Gal(Q/F)` は
   `Skeleton/IUTchI/InitialThetaData.lean` が記録している**オーバーバー脱落**（`Gal(F̄/F)` → `Gal(F/F)`）と
   同じ形をしており、**自明群と読めてしまう**。最初の目視対象はここにする。
2. **到達件数は pdftotext 由来**で、番号なしの依存を数え落とし、「cf.」型の案内を数え過ぎる。
   §2-0 で 1 件直したが、**同種の誤りが他にもありうる**。
3. **段 4 は実質「論文 2 本」**である（[Mzk1] の転写を含む）。
4. **Lemma 3.1 (iv) の [Serre] Ch IV §3.4 Lemma 3 を我々は読んでいない。** 「Frattini 型で自分で証明できる」は
   **見込みであって測定ではない**。段 1 に入った時点で最初に確かめること。

---

## 8. 実装状況

進捗の機械計数: **`node tools/genell-progress.mjs`**(`--json` で `ResearchPaper/genell-needed.json`)。

```
★ゴール進捗: [GenEll] の必要分 0 / 24 件 (0%)
  直接名指し 16 件 → 推移閉包 24 件 / GenEll 全 26 件
  ★条つき `.src` はあるが命題全体の `.src` が無いもの 1 件
     Lemma 3.1 — 条つき 3 個 ((i) / (ii) / (iii))
  §1 0/9   §2 0/1   §3 0/9   §4 0/5
```

★**分母は 24 件**(本書 §2 の「Cor 4.4 の到達 18 件」より大きい)。
需要側(IUTchI–IV・EtTh)が**直接名指しする 16 件**の閉包であり、
`Theorem 2.1` と §1 の 9 件を含むためである。両方の頂が分母に入っている。

### 段 0 —— 完了(2026-08-16)

| | |
|---|---|
| 目視 | 物理 p.13–14、200 dpi。`papers.json` に `GenEll` 登記(`notationRisk: high`) |
| 構造化 | `1_Structured/Arithmetic Elliptic Curves in General Position/`(index + lemma-3-1)。`check.mjs --structured` **PASS**(構造単位 112 件・逐語照合 112 件) |
| Lean | `lean/ABC3/Found/GenEll/Lemma31.lean`。**(i)(ii)(iii) を `sorry` 無しで証明**。`#print axioms` は標準 3 公理のみ |
| 逐語 | Lean コメント内の引用 6 件を PDF に対して照合、**NG 0**(すべて `layout` モードで一致) |
| 負の対照 | `lemma_3_1_ii_neg_control` —— `l = 2, 3` では全ての `λ` が `λ² = 1`。**`l ≥ 5` は飾りではない** |

★**`.src` は条つき 3 個のみ**。`Lemma 3.1` 全体の `.src` は (iv) が済むまで付けない
(器具の 2 値規則どおり。ゆえに進捗の数字は 0 のまま動かない——**これは正しい挙動**である)。

### ★実装して分かったこと(原文を読むだけでは出なかったもの)

1. **(i) の原文の道は `SL₂` の外に出る。** 原文は `γ ≝ (0 1 / 1 0)` を経由するが
   **det γ = −1** なので `γ ∉ SL₂(𝔽_l)`。原文はこれを注記していない。
   Lean では分解 `g = upper((a−1)/c)·lower(c)·upper((d−1)/c)` を明示する方が短い。
2. **(iii) に Sylow 論は要らない。** `W ≝ upper(−a/c)·M` の左上成分が 0 になることを使えば
   `W·upper(u)·W⁻¹ = lower(1)` が**行列の等式 1 本**で出る。
   原文の道(l-Sylow の個数 = l+1、正規化群 = Borel)は
   `|GL₂(𝔽_l)|` の計算と Borel の同定を要し、どちらも mathlib に無い。
3. **`l ≥ 5` が効くのは (ii) だけ**である(原文の証明も同じ)。
4. ★**器具の発見**: mathlib の `ring` は失敗すると `ring_nf` にフォールバックして
   **成功扱いになる**ため、`first | ring | …` は次の候補へ進まない。`ring1` を使う。
   本作業で 2 回これに当たった。

### ★段 1 の柱 1 を先に取った(同日)

段 0 の実装中に、(i) の証明が**体を一切使っていない**ことが分かった——
使うのは「`c` が可逆」か「`a` が可逆」の二択だけで、それは**局所環の定義そのもの**である。
そこで `mem_of_upper_lower_mem` を局所環 `[CommRing R] [IsLocalRing R]` へ一般化し、
系として次を得た(`sorry` 無し・標準 3 公理):

```lean
theorem closure_elementary_padic (l : ℕ) [Fact (Nat.Prime l)] :
    Subgroup.closure ((Set.range (upper : ℤ_[l] → SL(2, ℤ_[l]))) ∪
                      (Set.range (lower : ℤ_[l] → SL(2, ℤ_[l])))) = ⊤
```

★**`SL(2, ℤ_l)` は基本行列で代数的に生成される。位相を一切使っていない。**

さらに、その系として**還元の全射性**も (i) だけで出た:

```lean
theorem redPadic_surjective : Function.Surjective (redPadic l)
    -- redPadic l : SL(2, ℤ_[l]) →* SL(2, ZMod l) = SpecialLinearGroup.map PadicInt.toZMod
```

★`SL(2, 𝔽_l)` が `upper 1` と `lower 1` で生成され、その 2 つが明らかに持ち上がるから。
**Hensel の補題も位相も要らない。**

### ★★段 2 の入口で原典の記号の食い違いを発見した(2026-08-16)

§1 に入って最初に実装したのが `Definition 1.2, (ii)`(BD-class)である。
★**Arakelov 理論を要さない唯一の条**であり、しかも **abc の結論がこの語彙で書かれている**。

そこで**原典の記号 `≲` の向きが食い違っている**ことが判明した。**3 箇所すべて 260 dpi の PDF 目視**:

| 出典 | 目視した内容 |
|---|---|
| [GenEll] p.5, Def 1.2, (ii) | `α ≲ β` ⟺ ∃C, **`β(x) − α(x) ≤ C`**(= `β − α` が**上**に有界) |
| [GenEll] p.11, Thm 2.1, (i) | 表題が **“Effective Mordell/ABC/Vojta Conjecture”**、内容は `ht ≲ (1+ε)(log-diff + log-cond)` |
| [IUTchIV] p.54, Cor 2.3 | 同じ式に「i.e., the function `(1+ε)(…) − ht` **is bounded below** by a constant」<br>しかも `[cf. [GenEll], Definition 1.2, (ii)]` と明示的に引く |

★**ここまでが観測。以下は推論であり、区別して書く。**

**推論**: abc/Vojta の内容は `ht ≤ (1+ε)(…) + C`(= `(1+ε)(…) − ht` が**下**に有界)。
Thm 2.1 の表題と [IUTchIV] の言い換えは互いに整合し、abc の内容とも整合する。
Def 1.2, (ii) の印字だけがそれと逆になる。
→ **最も無理のない読みは「(a) で α と β が入れ替わっている(誤植)」**である。
★**これは推論であって観測ではない。**

**トリアージ**(`PLAN.md` §5-2): **① は読字については排除**(3 箇所とも PDF 目視、`pdftotext` 非経由)。
ただし「abc の内容は `ht ≤ (1+ε)(…) + C`」という**我々の数学的了解**が誤っている可能性までは
排除していない。② は該当しない。③ は上の推論が正しければ該当するが、
これは §5 が想定した「**飛躍**」ではなく**記号定義の誤植**であり、`Gap` の機構は当たらない。

★**独立確認: 取れた**(2026-08-16)。文脈を持たない子に、**結論を伏せたまま**
同じ 3 箇所を読ませた(260 dpi 全頁 ＋ 該当箇所 500 dpi 拡大)。子は各所を逐語で引いたうえで:

> **A4: 逆** —— 同じ差 `(1+ε)(log-diff + log-cond) − ht` を、
> (a) は「≤ C(上に有界)」と定義し、(c) は「bounded below by a constant(下に有界)」と
> 言い換えているから。

★これは `GAP` 型(否定的)の判定であり、`memory/challenger-audit-without-context.md` の
基準では `OK` より強い証拠である。**読字については独立に裏が取れた。**
残る共有前提は「abc の内容は `ht ≤ (1+ε)(…) + C` である」という数学的了解のみで、
そこは第三者の確認を経ていない。

### ★★Challenger 監査が実際に 1 件落とした —— `Definition 1.2, (ii)`(2026-08-16)

`BDClass.lean` にも同じ規律を当てた(依頼 3)。結果は **OK 4 / GAP 2 / UNVERIFIED 0**。
落ちたのは同一の 1 項目である:

> **GAP** —— 原文は「Finally, we observe that it makes sense to apply the notation
> “≳”, “≲”, “≈” **to BD-classes**.」と述べるが、Lean にあるのは関数レベルの両立性
> `bdle_congr` / `bdge_congr` の 2 本だけで、`BDClass F` 上の関係は 1 つも無く、
> **記号は BD-class に適用されていない**。

★**指摘は正しかった。** 親は「両立性(= well-defined 性の前提)を示したから適用できる」と
思っていたが、原文が言っているのは**適用そのもの**である。
`BDClass.le` / `BDClass.ge` / `BDClass.mk_eq_mk`(≈ の適用 = 類の等号)/
`BDClass.eq_of_le_of_ge` を足して埋め、**再監査(依頼 4)で `OK 1 / GAP 0`** を得た。

★**これが「完成宣言は文脈を持たない子に監査させる」規律の値打ちである**——
`memory/challenger-audit-without-context.md` が記録した `Prop 1.10` の実例と同型の
「親には構造上見えない欠落」が、`.src` を付けた直後の項目でもう 1 件出た。

★**運用上の学び**(子を使った 4 件で実測):
1. **報告手順を書かないと届かない**——依頼 1 では子が自力で `parent_send` を見つけたが、
   依頼 2 では黙って自分のセッション内に書いて終わった。**依頼文に報告手段を明記すること。**
2. **複数セッションが動いていると `wait_for_child` / `child_read` の出力が混ざる**
   (`[session switched → …]` が挟まる)。★**確実なのはファイル経由**——
   リポジトリ外の一時ファイルに書かせて親が読む方式に切り替えて解決した。

★**扱い**: どちらも型に出し、**両者が互いに逆であることを定理にした**——
`BDle`(印字どおり)/ `BDge`(abc の向き)/ `bdle_ne_bdge`(反例で別物と示す)。
**abc の主張を書くときは `BDge` を使う**という選択が、これで機械可読になった。

★**PLAN への含意**: `PLAN.md` §5 の類型(①モデル化 ②未構築 ③飛躍)に
**第 4 の類型「原典の誤植」が要る**。飛躍と違い falsifier は要らず、
**別の箇所での同じ記号の用法**が証拠になる。今回は同じ著者の別論文が証拠だった。

### ゲートの状態(2026-08-16)

`node tools/check.mjs` **PASS**(別セッションのビルドが止まった隙に 2 回実行):
`lake build` 成功 / `Found/` の `sorry` 0 / `.src` つき原典項目 **125 件**(locator 検証 114 件)。
1_Structured は構造単位 **114 件**・逐語照合 114 件で PASS。

★ゲートが実際に落とした例が 1 件ある: `definition-1-2.html` の初版で
HTML 実体参照(`&isin;` `&beta;` 等)を使ったところ、`check.mjs` の変換表に無いため
**S4 逐語照合が落ちた**(「先頭 27/72 文字まで一致」と場所まで印字された)。
リテラル文字に直して PASS。**器具は書き手の手癖を実際に捕まえる。**

## 8-2. ★★規模の結論(2026-08-16、実測)——24 件は 2 つの不在に律速される

§1 に入って測った結果、**残り 23 件はすべて、mathlib に無い 2 つの古典理論のどちらかに律速される**:

| 塊 | 件数 | 律速するもの | mathlib(探索範囲つき) |
|---|---:|---|---|
| §1(高さ) | 9 | **Arakelov 理論**——hermitian 計量つき算術直線束・`APic`・`ADiv`・`deg_F` | ★**0 件**(`Arakelov` / `arithmetic.*line bundle` で mathlib 全体を case-insensitive grep) |
| §2(Theorem 2.1) | 1 | §1 ＋ noncritical Belyi maps | Belyi **0 件**。ただし [Mzk1] は `0_Source` にある |
| §3(Galois 作用) | 9 | **l 捩れ点への Galois 表現**・Tate 曲線・半安定還元 | ★**0 件**(`EllipticCurve/*.lean` の `torsion` 8 件はすべて 2-torsion polynomial) |
| §4(素数の大きさ) | 5 | §3 | — |

★**あるもの**: `Mathlib/NumberTheory/Height/` **2006 行**(`mulHeight`/`logHeight`/Northcott/
`InfinitePlace` 経由の `archAbsVal`)。ただしこれは **ℙⁿ 上の Weil 高さ**であって、
**任意のスキーム上の算術直線束に付随する高さではない**。Def 1.1 の枠組みへは直接は載らない。

> ★**したがって B トラックは「IUT の転写」ではなく、
> Arakelov 理論と Galois 表現論という古典的数学を 2 本作る仕事である。**
> 規模は [FrdI] トラックと同程度かそれ以上——`PLAN.md` §9 の「見積りを Phase 2 の前に出さない」に従い、
> person-month/year の数字は出さない。**言えるのは「2 本の古典理論が丸ごと要る」までである。**

### ★この 2 本は**グラフに出ていなかった**——基礎理論の層を足した(2026-08-16)

`dependency-graph.html` を検索したところ、合流前は `Arakelov` **0 件** / `torsion` **0 件**だった。
既存グラフの 3 層(番号付き項目・§0 語彙・概念)のどれにも出ない——
★**概念層が拾えないのは、原文に定義文が無いから**である。原文は
「the Galois representation … associated to `E_L`」と**使う**だけで、定義文型では導入しない。

そこで `ResearchPaper/foundations.json` を作り、`tools/graph-layers.mjs` に
**基礎理論を出辺のない節点(= 真の葉、層 0)** として合流させた。9 節点・辺 27 本。
★張れなかった 3 件(`GenEll Definition 1.1` / `Corollary 4.3` / `Corollary 4.4`)は、
**本グラフの根が [IUTchIII] Corollary 3.12 だから**である——これらを直接引くのは IUTchIV。
**グラフに無いことは「要らない」ではない。**

★**この測定の価値**: 着手前は「§1 は mathlib に高さがあるから軽い」と見込んでいた(本書 §3-2)。
**実測はその見込みを外した**——高さはあるが**枠組みが違う**。
`PLAN.md` §8-1 の教訓(statement が言及する対象だけを見て証明の依存を見ない)と同型の誤りを、
今回は**着手前ではなく着手 1 日目に**捕まえたことになる。

## 9. 次の 1 手

**段 1 = `Lemma 3.1, (iv)` の残り(柱 2)。**

| 柱 | 内容 | 状態 |
|---|---|---|
| 1 | `SL(2, ℤ_l)` が `upper t` / `lower t` で生成される | ★**済**(上記) |
| 2 | 閉部分群 `J` が `upper t` / `lower t` を**全部**含むこと | **未** |

柱 2 が [Serre] Ch IV §3.4 Lemma 3 に当たる部分であり、pro-l 群の位相が要る。
**mathlib の実測(2026-08-16、探索範囲つき)**:

| 要るもの | mathlib | 探索範囲 |
|---|---|---|
| 還元 `SL₂(ℤ_l) → SL₂(𝔽_l)` | ★**ある**: `Matrix.SpecialLinearGroup.map (f : R →+* S)` | `LinearAlgebra/Matrix/SpecialLinearGroup.lean` |
| profinite 群 | **ある**(`ProfiniteGrp` 135 件) | mathlib 全体 grep |
| **pro-p 群** | ★**実質無い**(2 件) | mathlib 全体 grep(`IsProPGroup` / `pro-p`) |
| `ℤ_l` 上の合同部分群の族 | ★**無い**。`SL(2, ℤ)` 用のものが `NumberTheory/ModularForms/CongruenceSubgroups.lean` にあるだけ | 同上(`congruence.*subgroup`) |

→ 柱 2 は「壁の種類 1(mathlib の不在)」であり、**合同フィルトレーション
`Γ(l^n)` とその商 `≅ 𝔰𝔩₂(𝔽_l)` を自分で作る**ことになる。

### ★柱 2 を 4 段に割り、うち 1 段を実装した(2026-08-16)

(iv) の証明を解析して、[Serre] の補題の中身を次の 4 段に割った:

| 段 | 内容 | 状態 |
|---|---|---|
| 1a | `K_n / K_{n+1}` が**加法**群で、跡 0 の行列で書ける(代数部分) | ★**済**(`Found/GenEll/Sl2Congruence.lean`) |
| 1b-有限 | 合同核が `ℤ/l^{n+1}` の中で `p = l^n`、`p² = 0` を満たす | ★**済**(`pow_self_mul_self_eq_zero` / `congruence_mul` / `congruence_det`) |
| 1b-共役 | **共役作用が随伴作用に落ちる**(`g(1+p·A)g⁻¹ = 1+p·(gAg⁻¹)`)と、跡が保たれる | ★**済**(`conj_one_add_smul` / `trace_conj`。★`p² = 0` すら要らない) |
| 1b-位相 | `ℤ_l` の逆極限として `K_n` を定義する | **未**(位相が要る。★ここだけ) |
| 2 | `J ∩ K_n` の像が部分加群になる | **未** |
| 3 | **`𝔰𝔩₂(𝔽_l)` が既約**——部分加群は `0` か全体 | ★**済**(`Found/GenEll/Sl2Adjoint.lean`) |
| 4a | `Σ_{i<l} (1+D)^i = 0`(`D³ = 0`・標数 `l`・`l ≥ 5`) | ★**済**(`Sl2Congruence.lean`) |
| 4b | それを `x^l` の類の計算に接続する(`Ad(u) = 1 + D`、`D³ = 0` の同定) | **未** |

★**段 1a で取った代数の核**(`p² = 0` なる `p` に対して):

```lean
theorem mul_one_add_smul : (1 + p • A) * (1 + p • B) = 1 + p • (A + B)
theorem det_one_add_smul : (1 + p • A).det = 1 + p * (A 0 0 + A 1 1)
```

★**積が和に化けるのは `p²` の項が消えるからで、それ以外の理由は無い**。
そして 2 番目が「**なぜ跡 0 なのか**」を説明する——`det = 1` は `p · tr(A) = 0` と同値である。

★**段 4 の算術核は取った**(`Sl2Congruence.lean` の `WhereFiveEnters` 節):

```lean
theorem choose_one_two_three_eq_zero (l : ℕ) (hl : Nat.Prime l) (h5 : 5 ≤ l) :
    ((l.choose 1 : ℕ) : ZMod l) = 0 ∧ ((l.choose 2 : ℕ) : ZMod l) = 0
      ∧ ((l.choose 3 : ℕ) : ZMod l) = 0

theorem choose_three_three_ne_zero_mod_three :   -- ★負の対照
    ((Nat.choose 3 3 : ℕ) : ZMod 3) ≠ 0
```

★**負の対照が `l ≥ 5` の必要性を機械で固定する**——`l = 3` では `C(3,3) = 1 ≠ 0` なので
段 4 の消去が破れる。**`l ≥ 5` は飾りではない**。

★★**そして段 4 の代数は閉じた**:

```lean
theorem sum_range_choose_eq (n j : ℕ) :          -- hockey-stick(range 版)
    ∑ i ∈ Finset.range n, i.choose j = n.choose (j + 1)

theorem one_add_pow_of_cube_eq_zero (hD : D ^ 3 = 0) (i : ℕ) :
    (1 + D) ^ i = 1 + (i : R) * D + (i.choose 2 : R) * D ^ 2

theorem sum_range_one_add_pow_eq_zero (l : ℕ) (hl : Nat.Prime l) (h5 : 5 ≤ l)
    (hchar : (l : R) = 0) (hD : D ^ 3 = 0) :
    ∑ i ∈ Finset.range l, (1 + D) ^ i = 0        -- ★段 4 の核
```

★**`Σ_{i<l} (1+D)^i = 0`** ——これが「`x^l` の類が**持ち上げ方に依らない**」ことの中身である。
`l ≥ 5` は `C(l,3)` を消す 1 箇所だけで効く。

### ★★位相が要るのは最後の 1 歩だけ、と切り分けた(2026-08-16)

(iv) は `SL₂(ℤ_l)` の**閉**部分群についての主張だが、証明を 2 つに割ると:

| | 内容 | 位相 |
|---|---|---|
| (A) | `H ≤ SL₂(ℤ/l^n)` が `SL₂(𝔽_l)` へ全射なら `H = ⊤` | ★**要らない**(有限群論) |
| (B) | 閉部分群 `J ⊆ SL₂(ℤ_l)` は各 `mod l^n` の像で決まる | 要る(逆極限) |

★**(A) の帰納法に要る材料はすべて揃った**:

| 材料 | 宣言 |
|---|---|
| 合同核が加法群 | `mul_one_add_smul` / `slOfSmul_mul` |
| `det = 1 ⟺ 跡 0` | `det_one_add_smul` |
| `p = l^n` が `ℤ/l^{n+1}` で `p² = 0` | `pow_self_mul_self_eq_zero` |
| **共役作用 = 随伴作用**、跡を保つ | `conj_one_add_smul` / `trace_conj` |
| **核の元が `1 + l^n·A` と書ける** | `exists_matrix_of_castHom_eq_zero` |
| `𝔰𝔩₂(𝔽_l)` の既約性 | `sl2_adjoint_irreducible` |
| `Σ_{i<l}(1+D)^i = 0`(`l ≥ 5`) | `sum_range_one_add_pow_eq_zero` |

**残るのは (A) の帰納法を組むことと、(B) の逆極限**である。
★**位相を要するのは (B) だけ**——他はすべて有限段の代数で閉じている。

★**段 4 が `l ≥ 5` の本当の在り処である**。`x` を `upper 1` の持ち上げとすると
`x^l` の `K_1/K_2` における類は `E₁₂ + (Σ_{i<l} Ad(u^i))·A` で、
`u` が冪単なので `Σ_{i<l} Ad(u^i) = Σ_j C(l, j+1) N^j`(`N³ = 0`)。
`C(l,1) = l`、`C(l,2) = l(l−1)/2`、**`C(l,3) = l(l−1)(l−2)/6`** がすべて `≡ 0 (mod l)` に
なるのに **`l ≥ 5`** が要る(`l = 3` では `C(3,3) = 1 ≢ 0`)。
→ ★**`l ≥ 5` は「持ち上げ方に依らず `x^l` の類が `E₁₂` になる」ことの条件**である。
(この計算はまだ Lean にしていない。上の段 4。)

★**段 3 で分かったこと**: 既約性そのものに要るのは **`2 ≠ 0` だけ**であった。
`l ≥ 5` は要らない。仮定は `Submodule K` で、`AddSubgroup` で書くと
**スカラー倍で閉じないので `e` が取り出せない**(最初そう書いて失敗した)。

★**2026-08-16 に 2 経路を検討し、どちらも同じ壁に当たると測った**:

| 経路 | 要るもの | mathlib |
|---|---|---|
| (a) pro-l 群の Frattini 論法 | pro-p 群の理論 | **2 件**(実質不在) |
| (b) 合同フィルトレーション ＋ 随伴表現の既約性 ＋ 非分裂性 | 標数 l の modular representation theory | **1 件**(実質不在) |

`PLAN.md` §8-3 の stop-loss(「迂回路を 2 つ以上試して全部同じ壁」)に該当する。
★ただし (b) の中の「`𝔰𝔩₂(𝔽_l)` が `SL₂(𝔽_l)`-加群として既約(l ≥ 5)」は
**初等的な線形代数で書ける**ので、着手するならそこが入口である。

## 9-2. ★★24/24 への段階つき長期ゴール(2026-08-16 設定)

### ★まず前提: **スケルトン化では 24/24 にならない**

`tools/genell-progress.mjs` は **`Found/` の `.src` だけ**を数える
(`rel.split('/')[0] !== 'Found'` で除外。設計どおり——この数は「実装」を測る)。
★したがって `Skeleton/` を厚くしても分子は 1 も動かない。

**そのかわり副指標を足した**(同ファイル、2026-08-16):

```
★ゴール進捗: [GenEll] の必要分 0 / 24 件 (0%)
  ★副指標(**完了の数ではない**): statement が固定されている項目 1 件 / 手つかず 23 件
     Theorem 3.8 — Interface/GenEll/GaloisRep.lean / Skeleton/GenEll/GaloisImage.lean
```

★**動く数と動かない数を並べて見せる**ことで、「スケルトン化＝進捗」という誤解を封じる。

### ★段階(各段に**機械で判定できる**受理条件を付けた)

| 段 | 内容 | 受理条件(機械) | 2 理論が要るか |
|---|---|---|---|
| **S0** ✅ | 到達点を型で固定する | 副指標「statement が固定」が **24/24** | 要らない(語彙は `Interface` に posit) |
| **S1** ✅ | 2 理論なしで到達できる 3 件を実装 | `genell-progress` が **3/24** | ★**要らない** |
| **S2** | Arakelov の層 A・B・C を作り §1 を実装 | `genell-progress` が **12/24** かつ `ArithLineBundle` が `waiting` でない | Arakelov |
| **S3** | Galois 表現を作り §3 を実装 | `genell-progress` が **21/24** かつ `GaloisRep` が `waiting` でない | Galois 表現 |
| **S4** | §2・§4 を実装(Belyi と `Cor 4.4`) | `genell-progress` が **24/24** | 両方 ＋ Belyi |

★**S1 が「2 理論なしで到達できる上限」である**(下記)。S2 以降は理論を作らない限り 1 件も動かない。

### ★★S1 の中身 —— 2 理論なしで到達できるのは **3 件**(実測 2026-08-16)

| 項目 | なぜ 2 理論が要らないか | mathlib |
|---|---|---|
| `Lemma 3.1` | **純群論**。(i)(ii)(iii) は実装済み、残りは (iv) の pro-l 群 | Sylow・交換子・transvection あり |
| `Lemma 4.1` | ★**Chebyshev の評価 (i)(ii) は「仮説」として与えられている**——PNT を証明する必要が無い(`Remark 4.1.1` が「PNT の帰結」と明記して仮説の外に置く) | ★`Chebyshev.theta`(`θ(x) = Σ_{p≤x} log p`)が**ある** |
| `Lemma 4.2` | 「Some Elementary Estimates」。原文の証明が **2 行** | 初等 |

★**この 3 件は `Interface` を使わずに `Found/` へ載る**——だから空虚になりようがない。
★残る 21 件は `§1 → Arakelov` / `§3 → Galois 表現` / `§2・§4 → その上` で、
**1 件も 2 理論を迂回できない**(`ResearchPaper/foundations.json` の `neededBy`)。

### ★受理条件に「sorry 0 / axiom 0」だけを使わない理由

`PLAN.md` §1 が実演したとおり、**条件を取り違えた `structure` は
sorry 無し・公理依存ゼロのまま `0 = 1` まで証明できる**。
本セッションでも `Skeleton/GenEll/Heights.lean` の sorry を
`Interface` への輸入で 0 にした実例が出た。ゆえに各段の受理条件は

1. `genell-progress` の分子(= `Found/` の条なし `.src`)
2. `Found/` の `sorry` 0
3. `axiom` 0
4. ★対応する `Interface` が **`waiting` でない**
5. ★G2 の**非空虚 witness**

の **5 つを同時に**要求する。★4 と 5 が無いと 1〜3 は空虚になれる。

## 10. ★次にやるべきことの候補(どれも重い。ユーザーの判断が要る)

| | 内容 | 規模 | 効き方 |
|---|---|---|---|
| **A** | **Arakelov 理論の最下段**——`ADiv(F)` / `APrc(F)` / `deg_F`(GenEll p.3–4)。`NumberField.InfinitePlace`(70 件)と `HeightOneSpectrum` の上に載る | 中 | §1 の 9 件すべての土台。**Def 1.1 の完成には算術直線束が要るので、この段階では `.src` は付かない** |
| **B** | **`𝔰𝔩₂(𝔽_l)` の既約性**(l ≥ 5) | 小〜中 | 段 1 柱 2 (b) の入口 |
| **C** | **[Mzk1] Noncritical Belyi Maps の構造化**(9 ページ、`0_Source` にある) | 中 | Theorem 2.1 の外部依存を 1 つ潰す。Lean を触らない |
| ~~**D**~~ | ~~段 0 の Challenger 監査~~ | — | ★**実施した(下記)** |

### ★D の結果 —— Challenger 監査(2026-08-16、文脈を持たない子)

`memory/challenger-audit-without-context.md` の規律どおり、**親の判断を書かず**、
出力を `OK:` / `GAP:` / `UNVERIFIED:` の 3 値に固定して依頼した。

**判定: OK 3 / GAP 0 / UNVERIFIED 0。**

★**子は `.txt` を使わず PDF を 260 dpi で画像化して読んだ**(依頼文で「PDF が権威」とだけ伝え、
行列の順序が壊れることは**伝えていない**)。その上で
`upper 1 = !![1,1;0,1] = α`・`lower 1 = !![1,0;1,1] = β` と同定している——
★**`.txt` の罠(α と β が入れ替わる)を独立に回避した**ことになる。

★子が指摘した事実(判定には影響しないが記録する):
**(i) と (iii) は原文の `l ≥ 5` を落としている。**
これは仮定の削除(= 一般化)であり、`l ≥ 5` を代入すれば原文の主張がそのまま出る。
`PLAN.md` の A6 が危険とするのは**仮説の強化**であって、その鏡像である**弱化は安全**である。
実際 Lean は `5 ≤ l` 無しで型検査を通しており、**主張が真であることの証拠は証明そのもの**である。

★**なお `OK` は `GAP` より弱い証拠である**(同 memory)。そこで、より強い証拠が取れる
**否定的な問い**(記号 `≲` の向きが 3 箇所で整合するか)を、結論を伏せたまま別途ぶつけた。

## 9-3. ★★S0・S1 到達の記録と、S2 以降の壁の実測(2026-08-17)

### ★S0 到達 —— 24/24(statement を型で固定)

```
★副指標(**完了の数ではない**): statement が固定されている項目 24 件 / 手つかず 0 件
```

新設: `Interface` 4 本(`HeightTheory` / `AbcSetup` / `TateLocal` / `EllModuli`)、
`Skeleton` 4 本(`Section1`–`Section4`)、`Gap` 1 本(`BDDirection`)、
`1_Structured` 4 本(`section-1`–`section-4`.html、逐語照合 143 単位 PASS)。

★`Interface` 4 本は**公理を 1 つも持たない**(データと述語だけ)。
ゆえに `Prop 1.4 / 1.6 / 1.7`・`Lemma 3.2`・`Theorem 2.1` の結論は
**本当に `sorry` のまま残る**——`check.mjs` 冒頭 B5 の穴に落ちようがない。

### ★S1 到達 —— 3/24(`Lemma 3.6` / `Lemma 4.1` / `Lemma 4.2`)

★**S1 の中身は当初の見立てと 1 件違った。**
上の 9-2 の表は `Lemma 3.1` / `Lemma 4.1` / `Lemma 4.2` の 3 件としていたが、
**実測の結果 `Lemma 3.6`(An Elementary Estimate)が 4 件目**として加わり、
逆に **`Lemma 3.1` は (iv) が残るため到達しなかった**。

| 項目 | 状態 | 実測 |
|---|---|---|
| `Lemma 3.6` | ✅ 実装 | ★**純粋な実解析**。原文は極限で片づけるが、**極限では定数が出ない**ので `C₀ ≝ ((1+ϵ)/ϵ)^{1+ϵ}` を明示構成した |
| `Lemma 4.1` | ✅ 実装 | ★条件 (i)(ii) は**仮説**。原文末尾の WLOG は `θ` を `ℝ′_{>0}` から `ℝ_{>0}` へ延ばして**消した** |
| `Lemma 4.2` | ✅ 実装 | ★`H+1 ≤ 2^H` から原文より強い中間評価を経由(主張は原文どおり `3h/2`) |
| `Lemma 3.1` | ❌ **(iv) が残る** | (i)(ii)(iii) は実装済み。(iv) は [Serre] IV, §3.4, Lemma 3(`0_Source` に無い) |
| `Remark 4.1.1` | ❌ | ★後半が**素数定理そのもの**。mathlib に無い(`Chebyshev.lean` の全定理名を確認) |

★**2 理論なしで到達できる上限は 4 件**(3 件ではない)。現在 3 件。

### ★★S2 以降の壁 —— 何が足りないかを件数で書く

★S2(12/24)以降は、**mathlib に 0 件の理論を作ること**が受理条件に入っている。

| 段 | 追加で要る理論 | mathlib 実測(2026-08-16/17) | 公開プロジェクト |
|---|---|---|---|
| **S2** | 算術直線束 `APic(X)`、複素解析空間 `X^arc`、Faltings 高さ | `Arakelov` **0 件** / `arithmetic line bundle` **0 件** / `LineBundle` **0 件** / `analytification` **0 件** / `GAGA` **0 件** / `complex analytic space` **0 件** | 無し(実測) |
| **S3** | l 進 Tate 加群 `M_l(E)` と Galois 作用、Tate 曲線、半安定還元 | `E[n] ≅ (ℤ/n)²` **無し** / Galois 作用 **無し** / Tate 曲線 **無し** | ★**FLT プロジェクト**(Buzzard、EPSRC EP/Y022904/1、2029-09 まで)ブループリント §2.5・§3 と重なる |
| **S4** | Belyi 写像(noncritical) | `Belyi` **0 件** | 無し(実測) |

★★**したがって S2・S3・S4 は「作業を続ければ届く」段ではない。**
それぞれが**独立した基礎理論の構築**であり、S3 に至っては
**国家助成つきの多年度プロジェクトが 2029 年を目標にしている範囲**と重なる。

★**この文書は S2–S4 を諦めるためのものではない。**
受理条件を機械判定にしてあるのは、**部分的な前進を数えられるようにする**ためである:

- `ArithLineBundle` / `GaloisRep` の `waiting` が**何を待っているか**は型に書いてある
- 各段の分子は `node tools/genell-progress.mjs` がいつでも出す
- ★**「まだ 3/24」と「もう 3/24」の両方が同じ数で言える**のがこの指標の役目である

### ★次の 1 手(2 理論を待たずにできる唯一のもの)

**`Lemma 3.1, (iv)`** —— これだけが 2 理論なしで 4/24 へ進む道である。

★段 (A)(代数)の部品は**既に実装済み**:

| 部品 | ファイル | 役割 |
|---|---|---|
| `congruence_mul` | `Sl2Congruence.lean` | `K_n/K_{n+1}` が**加法**群 |
| `congruence_det` | 同 | 行列式が跡を測る |
| `conj_one_add_smul` | 同 | 共役作用 = **随伴**作用 |
| `trace_conj` | 同 | 随伴は跡を保つ(`𝔰𝔩₂` が不変) |
| `exists_matrix_of_castHom_eq_zero` | 同 | 合同核の元は `1 + l^n·A` の形 |
| `sl2_adjoint_irreducible` | `Sl2Adjoint.lean` | `𝔰𝔩₂(𝔽_l)` の**既約性**(`2 ≠ 0` だけで足りる) |
| `sum_range_one_add_pow_eq_zero` | `Sl2Congruence.lean` | `x^l` の計算(`M = 0` の場合を潰す段) |

★残るのは (a) 有限段の帰納法をこれらで組むこと、(b) `ℤ_l` の位相で閉包を取ること。
**(b) が「位相が要る最後の 1 歩」**である(§9 の切り分けを参照)。

## 9-4. ★★`Lemma 3.1` 完成の記録と、S2 の分解(2026-08-17)

### ★4/24 —— 2 理論なしで届く上限に到達した

```
★ゴール進捗: [GenEll] の必要分 4 / 24 件 (17%)
  §3   2/ 9  ██·······     §4   2/ 5  ██···
```

| 項目 | 実装 | 実測して分かったこと |
|---|---|---|
| `Lemma 3.1` | `Lemma31.lean` ＋ `Sl2Level.lean` ＋ `Sl2Padic.lean` | ★(iv) の [Serre] Ch IV §3.4 Lemma 3 は `0_Source` に無く mathlib にも無い——**自分で証明した**(1094 行) |
| `Lemma 3.6` | `Elementary.lean` | ★原文は極限で片づけるが**極限では定数が出ない**。`C₀ ≝ ((1+ϵ)/ϵ)^{1+ϵ}` を明示構成 |
| `Lemma 4.1` | `PrimesOfSize.lean` | ★原文末尾の WLOG を、`θ` を `ℝ′_{>0}` から `ℝ_{>0}` へ延ばして**消した** |
| `Lemma 4.2` | `Elementary.lean` | ★`H+1 ≤ 2^H` で原文より強い中間評価(主張は原文どおり) |

### ★副次的に得たもの(いずれも mathlib に無い)

| 宣言 | 内容 |
|---|---|
| `sl2_commutator_eq_top` | ★**`SL₂(𝔽_l)` は完全群**(`l ≥ 5`)。`layer0-triage.md` が「群の完全性は mathlib に無い」と記録していたもの |
| `subgroup_eq_top_of_redPow_surj` | `SL₂(ℤ/l^n)` の合同部分群論(**位相不要**、813 行) |
| `sl2_padic_eq_top_of_isClosed` | [Serre] Ch IV §3.4 Lemma 3 の `SL₂` 版 |
| `sum_adU_eq_zero` | `Σ_{i<l^n} Ad(u^i) = 0` on `𝔰𝔩₂(𝔽_l)`——★演算子環(非可換)を避けて**成分で**取った |

### ★★残り 20 件は 1 件も 2 理論を迂回できない(再確認)

| 節 | 残り | 律速 |
|---|---:|---|
| §1 | 9 | ★**Arakelov**——`APic(X)`(算術直線束)と `X^arc`(複素解析空間) |
| §2 | 1 | Arakelov ＋ **Belyi** |
| §3 | 7 | **Tate 曲線**・**Faltings 高さ**(= Arakelov) |
| §4 | 3 | 上の全部(`Theorem 3.8` 経由) |

★例外は `Remark 4.1.1` の 1 件で、これだけは別種の壁——**素数定理**である
(mathlib の `Chebyshev.lean` は Chebyshev 型の評価を持つが `θ(x)/x → 1` は無い)。
公開プロジェクト `PrimeNumberTheoremAnd` が持つが、
★**使う前に clone して `sorry` を数えること**(`lean-ecosystem.json` の規律)。まだ測っていない。

### ★★S2(12/24)に何が要るか —— 分解して測った

S2 の受理条件は「12/24 かつ `ArithLineBundle` が `waiting` でない」。
`Interface/GenEll/HeightTheory.lean` の `waiting` が名指ししているものを分解する:

| # | 作るもの | mathlib | 依存 |
|---|---|---|---|
| A1 | **差積イデアル → 算術因子** `δ_F ∈ ADiv(F)` | ★`IsDedekindDomain.differentIdeal` が**ある** | `ADiv`(実装済み) |
| A2 | Cartier 因子の `Spec(O_F)` への引き戻し | scheme 論はあるが算術因子との橋は無い | scheme |
| A3 | **複素解析空間** `X^arc` とエルミート計量 | ★**0 件** | 複素解析 ＋ GAGA |
| A4 | **算術直線束** `APic(X)` | ★**0 件** | A3 |
| A5 | 高さ `ht_L̄` と `Proposition 1.4` | ★**0 件** | A4 |

★★**A1 だけが今すぐ取れる**——mathlib の `differentIdeal` と我々の `ADiv` が既にあるからである。
★A3 は**複素解析空間の理論そのもの**であり、mathlib に着手の形跡が無い(2026-08-16 実測)。
**S2 はこの 1 点で律速される。**

### ★この文書の役目(再掲)

★**「まだ 4/24」と「もう 4/24」の両方が同じ数で言える**ようにしてある。
`node tools/genell-progress.mjs` はいつでもその数を出し、
`node tools/check.mjs` は `Interface` の `waiting` が何を待っているかを並べる。

## 9-5. ★★公開プロジェクトの実測(2026-08-17)—— §9-3 の記述を 1 つ訂正する

規律どおり **clone して `sorry` を数えた**(`ResearchPaper/lean-ecosystem.json` に全文)。

### ★訂正: 「S3 は FLT と重なる」は**予定の話**であって実装の話ではなかった

§9-3 の表に「S3 … 公開プロジェクト: ★**FLT プロジェクト**…ブループリント §2.5・§3 と重なる」と書いた。
★**実測すると、我々が要る土台は FLT にも無い。**

| ファイル | 実測(2026-08-17) |
|---|---|
| `FLT/EllipticCurve/Torsion.lean` | ★**124 行中 sorry 10 件**。`n_torsion_card` / `n_torsion_dimension`(= `E[n] ≅ (ℤ/n)²`)が**いずれも sorry** |
| `FLT/TateCurve/TateCurve.lean` | ★**20 行**。Tate 曲線の理論ではなく入口の宣言だけ |
| `FLT/GaloisRepresentation/` | 611 行・sorry 9 件。★**Frey 曲線の 3 進表現に特化**。一般の `Gal(ℚ̄/L) → GL₂(ℤ_l)` ではない |

★toolchain も `v4.34.0-rc1`(我々は `v4.31.0-rc2`)で差が大きい。
★★**S3 は「FLT を待てば済む」段ではない。**

### ★★`Remark 4.1.1` は埋まりうる —— ただし依存の判断はユーザーのもの

`PrimeNumberTheoremAnd/Consequences.lean:177` に
**`chebyshev_asymptotic : θ ~[atTop] id`** がある。

★★**`θ` は mathlib の `Chebyshev.theta` そのもの**である
——我々が `Found/GenEll/PrimesOfSize.lean` で使っているものと**同一**で、変換が要らない。

★推移的 import 閉包は**プロジェクト内 9 ファイルだけ**、その sorry は**合計 2 件**で、
どちらも我々の経路の外にある(`prelim_decay_2` は**どこからも使われず**、
`prelim_decay_3` を使う `decay_alt` を使うのは閉包外の `IEANTN/CH2` だけ)。

★**ただし決定的ではない**——grep による依存追跡であって `#print axioms` ではない。
確定には build が要るが、toolchain が違う(`v4.32.2`)。

★★**依存に加えてはいない。** build 時間・保守・toolchain の巻き添えが伴うので、
これは**ユーザーの判断**である。採用すれば `genell-progress` は 4/24 → **5/24**。

### ★この測定が変えたこと

| | 測定前 | 測定後 |
|---|---|---|
| S3 の見通し | 「FLT と重なるので独立に作ると重複投資」 | ★**FLT にも無い。重複投資の心配より『誰も作っていない』方が問題** |
| `Remark 4.1.1` | 「素数定理が要るので不可能」 | ★**外部依存を 1 つ足せば可能**(判断はユーザー) |
| 到達可能な上限 | 4/24 | 4/24(依存を足せば **5/24**) |

★**S2 の律速(複素解析空間 `X^arc`)はどちらのプロジェクトにも無い。** ここは変わらない。

## 9-6. ★★基礎理論側の到達点(2026-08-17)—— 何が動き、何が動かなかったか

ユーザーの指示は「**基礎理論の構築を続ける**」。以下はその区間の記録である。

### ★★主指標は動いていない —— 4/24 のまま

★**これが最も重要な事実である。** 以下に並べるものはすべて `sorry` 0・公理 0 で、
ゲートも PASS だが、**24 件のうち 1 件も新たに完成していない**。
`Definition 1.5` が (ii)(iii) だけ済んで (i)(iv) で止まる、という形が典型である。

### ★過剰申告を 1 件出して、自分で直した

`deg_adivRed_le.src` に**条なしの `"Proposition 1.6"`** を付け、
`genell-progress` が **4/24 → 5/24 に誤って動いた**。
★`deg_adivRed_le` は `Proposition 1.6` の**非アルキメデス側だけ**である。

★★**機械は止めてくれなかった**——`check.mjs` は `sectionId` の実在しか見ず、
「その定理が本当にその項目の全体か」は見られない(冒頭 A 群と同じ性質)。
★**数が動いたことに気づいて内訳を追ったのが唯一の検出経路**だった。

> ★教訓: **数が動いたら、なぜ動いたかを必ず内訳で確かめる。**

### ★取ったもの(すべて `Found/`、`sorry` 0・公理 0)

| ファイル | 内容 | 原文の対応 |
|---|---|---|
| `ArithDivHom.lean` | `principalADiv` が準同型、`APrc(F)` は**像そのもの** | p.4 地の文 |
| `ProductFormula.lean` | ★★**積公式** → `deg` が `APic(Spec(O_F)) → ℝ` に降りる | p.4「determines a homomorphism APic(Spec(OF)) → R」 |
| `LogDiff.lean` | 差積イデアル → 有効算術因子 | `Definition 1.5, (iii)` |
| `LogDiffValue.lean` | ★★**`log-diff = log\|disc F\|/[F:ℚ]`**(値まで) | 同 |
| 同 | 導手の評価のイデアル側 | `Proposition 1.6` の非アルキメデス側 |
| `Conductor.lean` | `(−)_red` の冪等性・有効性・**次数不等式** | `Definition 1.5, (ii)` / `Proposition 1.6` |
| `NorthcottRat.lean` | ★**`ℚ` 上の Northcott 性** | `Proposition 1.4, (iv)` の基底 |

### ★★符号のバグを 1 件見つけた —— 「sorry 0 では足りない」の実例

`deg(ADiv(f)) = 0` を取りにいって符号が合わず、
**mathlib の `intValuationDef` が `exp(−count)` を返す**ことが判明。
私の `ordv` は原文の `ord_v` の**符号違い**だった。

★★**当該ファイルは最初から `sorry` 0・公理 0 で、それでも中身が誤っていた。**
`deg` を使う定理を 1 つも証明していなかったので、**型検査は通り続けていた**。
★`ordv_algebraMap_eq_count`(整元では `ord_v` = 重複度)を負の対照として足した。

### ★★mathlib の実測 —— 「分野がある」と「要る定理がある」は別

| 我々が要るもの | mathlib | 備考 |
|---|---|---|
| 高さ(`mulHeight` / `logHeight`) | ★**ある**(5 ファイル 2006 行) | |
| **Northcott 性**(= `Prop 1.4, (iv)`) | ★★**無い**。`Northcott.lean` が **TODO** と書き、基底インスタンスが 0 件 | ★`ℚ` の場合を本セッションで取った |
| 積公式 | ★**ある**(`prod_abs_eq_one`) | そのまま使えた |
| 差積イデアルのノルム = 判別式 | ★**ある**(`absNorm_differentIdeal`) | そのまま使えた |
| **無限素点の相対次数和** `Σ_{w\|v}[K_w:F_v]=[K:F]` | ★★**無い**(絶対版 `sum_mult_eq` のみ) | ★`degNormalized` の底変換不変性のコストを押し上げる |
| 根基の素因子分解 `radical I = ∏_{v∣I} v` | 無い | ★★**要らなかった**(下記) |

### ★「正面から解く必要は無かった」が 2 件増えた(累計 5 件)

4. **`Lemma 3.1, (iv)` に位相的閉包は要らなかった** —— `J.comap toGL` が既に閉。
5. **導手の評価に根基の素因子分解は要らなかった** —— `I ≤ radical I` だけで足りた
   (200〜300 行を見込んでいたが 10 行)。

★既存 3 件: round-trip 補題 / `charP_iff_prime_eq_zero` の探索範囲 / ノルムのダイヤモンド。

### ★★§1 は一枚岩ではない(4 例目まで来た)

| 箇所 | `X` が要るか |
|---|---|
| `Definition 1.5, (iii)`(log-diff) | ★**要らない**——`δ_x` は最小定義体 `F` だけの関数 |
| `deg` の商への降下 | ★`X = Spec(O_F)` で閉じる |
| `Proposition 1.6` の**非アルキメデス側** | ★**要らない**——算術因子で閉じる |
| `Proposition 1.4, (iv)` の**`ℚ` の場合** | ★**要らない**——有理数の高さで閉じる |
| 上記以外(`Definition 1.1`・`Example 1.3, (ii)`・`Prop 1.6` のアルキメデス側 …) | ★**要る**(`X^arc`) |

★★**したがって「§1 は Arakelov 待ち」は粗すぎる。**
**部品の過半は取れており、止まるのは `X^arc` と scheme の 2 点である。**

### ★次に取れるもの(と、その障害)

| 候補 | 障害 |
|---|---|
| `degNormalized` の底変換不変性 | ★**無限素点の相対次数和が mathlib に無い**(上記)。自前で要る |
| `Proposition 1.7` の証明の核(局所体の分岐＋Kummer) | ★局所体の差積イデアルの整備。`ClassFieldTheory` は該当箇所に sorry 11 件(2026-08-14 実測) |
| `Definition 1.5, (iv)` の引き戻し | scheme 論 |
| `Proposition 1.4, (iv)` の一般次数版 | 数え上げ(`ℚ` の場合は済) |

## §9-7 2026-08-17 —— `§1` の底を 4 箇所測り直した(3 箇所で自分の判断が誤っていた)

★本日の作業は「基礎理論の構築を続ける」というユーザの選択に従ったものである。
**`genell-progress` は 4/24 のまま動いていない**——それは想定どおりで、
本節が記録するのは**指標ではなく地図の訂正**である。

### 取ったもの(すべて `Found/`、sorry 0、`#print axioms` は標準 3 つのみ)

| ファイル | 中身 |
|---|---|
| `FinitePlaceRel.lean` | `hosComap`(素点の引き戻し)/ `hosFiber` / **`Σ_{w\|v} e·f = [K:F]`** |
| `BaseChange.lean` | `baseChange : ADiv F → ADiv K` / **`degNormalized` の底変換不変性** |
| `LogDiffTower.lean` | **`log-diff(F) ≤ log-diff(K)`**(= 原文が "minimal" と書く理由) |
| `CartierPullback.lean` | `IsEffectiveCartier` / `pullbackIdeal` / `conductorADiv` / **`log-cond`** |
| `MinField.lean` | **`F_min`** と**極小性**(`minField_le_of_factors`) |

### ★★自分の判断が誤っていた 3 箇所

1. **「`Cartier` が mathlib に 0 件 ⟹ 底が抜けている」は誤り。**
   名前は確かに 0 件だが、中身(因子の引き戻し)は
   `IdealSheafData.comap` + `equivOfIsAffine` + `ΓSpecIso` で**そのまま届いた**。
   ★**名前で測って「無い」と結論したのが誤りの型である。**
   同日中に `blocked-leaves.json` を書き換えて取り消した。

2. **`Definition 1.5, (i)` に scheme 論的像の理論は要らなかった。**
   `Scheme.SpecToEquivOfField : (Spec K ⟶ X) ≃ Σ x, (κ(x) ⟶ K)` が
   mathlib に**ある**ので、構成と分解が同時に手に入り、
   極小性は**単射性を 1 回**使うだけで出た。

3. **`degNormalized` の不変性は「`Interface` の仮説」ではなく定理にできた。**
   `ArithDiv.lean` は「不変性そのものは本ファイルでは示していない」と
   書いたままだった。有限側(基本等式)と無限側(相対次数和)を
   別々に組んで、`ADiv` の層では定理になった。
   ★★**ただし `X` 上の直線束の高さの不変性は別物**で、
   `Interface/GenEll/ArithLineBundle.lean` の `base_change_invariant` は
   **消えていない**。混同しない。

### ★★★止まったのは、いちばん素朴に見えた条だった

`Definition 1.5` は (i)(ii)(iii)(iv) の 4 条。本日 (i)(iii)(iv) が揃った。
**残ったのは (ii)** ——「正規ネーター scheme の正則部分に含まれる有効 Cartier 因子 `E`
について `E_red` もまた有効 Cartier 因子である」。

★中身は **Auslander–Buchsbaum(正則局所環は UFD)** に帰着する。
★★**mathlib 実測: 0 件。** `IsRegularLocalRing` は `RegularLocalRing/Defs.lean`
**1 ファイルだけ**で、定義と基本補題しかない。

★★**したがって `Definition 1.5` には条なしの `.src` を付けない。**
付ければ `genell-progress` は 4→5 に動くが、それは過剰主張である。
(本日、同種の過剰主張を `deg_adivRed_le` で 1 度やって取り消している。)

### ★「正面から解く必要は無かった」が 2 件増えた(累計 7 件)

6. **因子の引き戻しに Cartier 因子の理論は要らなかった** —— `IdealSheafData.comap` で足りた。
7. **`F_min` に scheme 論的像は要らなかった** —— `SpecToEquivOfField` で足りた。

### ★§1 の「一枚岩ではない」表を更新

| 箇所 | 状態(2026-08-17 夜) |
|---|---|
| `Definition 1.5, (i)` minimal field | ★**実装済み**(極小性つき) |
| `Definition 1.5, (ii)` `E_red` が Cartier | ★**Auslander–Buchsbaum 待ち** |
| `Definition 1.5, (iii)` log-diff | ★**実装済み**(値・tower・極小性の理由まで) |
| `Definition 1.5, (iv)` log-cond | ★**実装済み**(負の対照 `logCond_top` つき) |
| `degNormalized` の底変換不変性 | ★**実装済み** |
| `Proposition 1.6` 非アルキメデス側 | ★**実装済み**(導手 ≤ 引き戻した因子の次数) |
| `Proposition 1.6` の右辺(**高さ**) | ★`Definition 1.1` → `X^arc` 待ち |
| `Definition 1.1` 算術直線束 | ★`X^arc`(複素解析空間)待ち |

★★**§1 で残っている障害は 2 つだけになった: `X^arc` と Auslander–Buchsbaum。**

### §9-7 追記(同日、後半)—— さらに 4 本

| ファイル | 中身 |
|---|---|
| `PrimesOfSize.lean` (追記) | `exists_xeps_cond_ii` —— `Remark 4.1.1` の**第 1 項** |
| `RadicalPrincipal.lean` | ★**mathlib の TODO**: `Ideal.radical (span {a}) = span {radical a}` |
| `LogDiffFinite.lean` | ★**log-diff の Northcott 性**(次数を切れば体は有限個) |
| `HeightADiv.lean` | ★★★**`logHeight₁(x) = deg_F(polarADiv x)`** |
| `NorthcottNF.lean` | 数体の**整元**についての Northcott 性 |

#### ★★「高さ = 算術因子の次数」を `X^arc` 抜きで取った

`HeightADiv.lean` の `logHeight₁_eq_deg_polarADiv` は、
原文の「高さ = 引き戻した算術直線束の次数」という見方の
**最も単純な場合**(`X = ℙ¹`、`L = O(1)`)を、
**`X^arc` を一切使わずに**取ったものである。

★★**一般の `X` には届かない**——一般には `L` の計量が `X^arc` 上の
エルミート計量として与えられねばならない。**混同しない。**

★副産物: `ordv` の符号(本日一度直した)の**独立な検算**になった。
符号が逆なら極因子ではなく零因子になり、高さが負になる。

#### ★半分だけ取れたものを 3 つ、条つきで記録した

| 項目 | 取れた側 | 残った側 |
|---|---|---|
| `Remark 4.1.1` | 条件 (ii) の初等性 | **素数定理**(`θ(x)/x → 1`) |
| `Definition 1.5, (ii)` | UFD での主イデアルの根基 | **Auslander–Buchsbaum** + scheme 化 |
| `Proposition 1.4, (iv)` | 数体の**整元**の Northcott 性 | **分母イデアルのノルム評価** |

★★いずれも `.src` を**条つき**にしてあり、`genell-progress` は **4/24** のままである。
**半分取れたことを全部取れたと読まない。**

#### ★残った 3 つの穴の性質は違う

- **素数定理**: 公開プロジェクト `PrimeNumberTheoremAnd` に**ある**(sorry 無しに見える)。
  依存を足すかどうかの**判断**の問題であって、数学の問題ではない。
- **Auslander–Buchsbaum**: mathlib に無く、ホモロジー代数を要する。**数学の作業**。
- **分母イデアルのノルム評価**: 道筋は明確
  (`N(𝔡_x) = ∏ᶠ_v max(|x|_v,1) ≤ H(x)` と `Ideal.absNorm I ∈ I`)。
  **手間**の問題。★`HeightADiv.lean` の `finsum_posLog_finitePlace` が既にその半分である。

### §9-7 追記(同日、夜)—— ★★★数体上の Northcott 性が**完成した**

昼に「残った 3 つの穴」の 1 つとして
「`Proposition 1.4, (iv)` の**分母イデアルのノルム評価**(**手間**の問題)」
と書いた。★**その日のうちに埋めた。**

| ファイル | 中身 |
|---|---|
| `OrdvIntegral.lean` | `v(f) ≤ 1 ⟺ ord_v(f) ≥ 0`、`ord_v ≥ 0` なら整元 |
| `DenominatorBound.lean` | `denIdeal` / `N(𝔡_x) = ∏ᶠ_v max(\|x\|_v,1)` / **`{y : K \| H(y) ≤ B}` は有限** |

★★**mathlib の穴を 1 つ塞いだ**——`Height/Northcott.lean` は
条件つき instance を 1 つ持つだけで **`Northcott (mulHeight₁ (K := K))` の
基底 instance を 1 つも持たない**。それを与えた。
これが入ると mathlib 側の条件つき instance が発火して
`Northcott (logHeight₁ (K := K))` も自動で付く。

#### ★機構——高さの有限素点側が、そのまま分母イデアルのノルムだった

`HeightADiv.lean` で `max(|x|_v,1) = q_v^{max(−ord_v(x),0)}` を取っていたので、
**分母イデアル `𝔡_x ≝ ∏_v v^{max(−ord_v(x),0)}` のノルムが高さの有限素点側そのもの**
になった。`Ideal.absNorm I ∈ I`(mathlib)で、分母として使える
**具体的な自然数** `d ≝ N(𝔡_x) ≤ H(x)` が手に入る。

★★**高さの「積の形」(`H(d·y) ≤ H(d)·H(y)`)は一切要らなかった**——
共役の絶対値を直接評価すれば済んだ(`‖φ(d·y)‖ = d·‖φ(y)‖`)。
★「正面から要ると思ったものが要らなかった」の **8 例目**。

#### ★残った穴は 2 つになった

| 穴 | 性質 |
|---|---|
| **素数定理**(`Remark 4.1.1` の第 2 項) | 公開プロジェクトに**ある**。**判断**の問題 |
| **Auslander–Buchsbaum**(`Definition 1.5, (ii)`) | mathlib に無い。**数学の作業** |

★そして `§1` 全体を止めている本命は依然として **`X^arc`**(複素解析空間)である。

### §9-7 追記(同日、深夜)—— ★★★高さの拡大公式と、mathlib の TODO 3 件

| ファイル | 中身 |
|---|---|
| `NorthcottProj.lean` | ★**mathlib の TODO 2 件**: `Northcott (Projectivization.mulHeight/logHeight)` |
| `HeightExtension.lean` | ★★★**`H_L(x) = H_K(x)^{[L:K]}`** と **絶対高さ** `h(x) ≝ logHeight₁/[K:ℚ]` |

#### ★★★今朝の道具が夜に効いた

`HeightExtension.lean` の証明は、**今朝 `degNormalized` の底変換不変性のために作った
2 つの相対和をそのまま使う**:

- `sum_mult_comap_eq`(`Σ_{w|v} mult_w = mult_v·[L:K]`、`InfinitePlaceRel.lean`)
- `sum_ramification_inertia_hos`(`Σ_{w|v} e·f = [L:K]`、`FinitePlaceRel.lean`)

★★どちらも「**素点の族を `[L:K]` 倍に数える**」という同じ形だからである。
**朝の段階でこの再利用は見えていなかった。**

★`hosComap_liesOver` を instance にしておいたことも効いた——
mathlib の `HeightOneSpectrum.valuation_liesOver`(`ord_w(x) = e·ord_v(x)`)が
インスタンス引数でそれを要求する。

#### ★★`e` と `f` が両方要る構造が 2 度現れた

| 場所 | `e` の出所 | `f` の出所 |
|---|---|---|
| `BaseChange.lean`(算術因子) | 底変換の係数 | 重み `log q_w` |
| `HeightExtension.lean`(高さ) | `ord` の関係 | 剰余体の大きさ |

★同じ構造が**別の量について 2 度**現れた。偶然ではなく、
どちらも「素点を数える」量だからである。

#### ★本日 mathlib に無かったものを 4 つ与えた

1. `Northcott (mulHeight₁ (K := K))` の**基底 instance**(数体一般)
2. `Northcott (Projectivization.mulHeight)`(mathlib が TODO と明記)
3. `Northcott (Projectivization.logHeight)`(同上)
4. `Ideal.radical (span {a}) = span {radical a}`(mathlib が TODO と明記)

★および `HeightOneSpectrum.comap`、無限素点の相対次数和、
高さの拡大公式(いずれも検索して 0 件だったもの)。

#### ★★それでも `genell-progress` は 4/24 である

★**高さ理論の側は厚くなったが、`ht_L̄`(算術直線束の高さ)には
依然として `X^arc` が要る。** `Definition 1.1` が空いている限り
`§1` の項目は埋まらない。**混同しない。**

## §9-8 2026-08-17 —— [NCBelyi] 線を開けた(`Theorem 2.1` の枝)

★[GenEll] `Theorem 2.1` は `[NCBelyi] Theorem 2.5`(非臨界 Belyi 写像の存在)を引く。
その計算の中心が `Lemma 2.1`、受け皿が `Lemma 2.2` である。

### ★★260 dpi 目視で pdftotext の欠落を 2 つ発見した

物理 p.2(`Lemma 2.1`)を目視して:

1. 条件 (iii) の `α ≠ 0, r, 1, ∞` の **`≠` が `=` に化ける**
2. **`f′(x)` のプライム `′` が落ちる**(`f (x) = 0` に見える)。条件 (d) の `≠` も同様

★★**目視せずに pdftotext だけで読んでいたら、
条件 (iii) を「`α = 0, r, 1, ∞` なら `α > 1`」と真逆に転写していた。**
`papers.json` の NCBelyi に登記し、`verifiedPages` に p.2 を追加した。
★これで NCBelyi の既知の欠落は 4 件(`ℚ̄` の上線 / `∩` / `≠` / `′`)。

### 実装した部品(すべて sorry 0、`#print axioms` は標準 3 つのみ)

| 部品 | ファイル |
|---|---|
| `Lemma 2.1` (b) 臨界点・(c) `f(β) ∉ f(S)`・(d) の 4 本 | `Found/NCBelyi/Lemma21.lean` |
| λ 正規化(`λ ≝ 1/α₂`)・`r = m/(m+n)`・`\|S′\| < \|S\|`・基底段 | `Found/NCBelyi/Normalize.lean` |
| 相対版の合成(帰納法が実際に使う形) | `Found/NCBelyi/BelyiComp.lean` |

★残るのは**帰納法そのものを回す段**だけである。

### ★★原文の粗描から特定したもの

`Lemma 2.2` の証明は原文で 6 行である。そこから明示できたもの:

| 原文の書き方 | 特定した中身 |
|---|---|
| 「for some appropriate positive rational number λ」 | **`λ ≝ 1/α₂`**(2 番目に小さい元) |
| 「so long as `\|S\| ≥ 4`」 | **`α₂` の存在条件**そのもの |
| 「`\|S′\| < \|S\|`」 | 根拠は **`f(0) = f(1) = 0`** |
| 「composing the resulting morphisms」 | Belyi 型どうしではなく**相対版**(★下記) |
| (書かれていない)`\|S\| ≤ 3` | **1 次式 `x ↦ c·x`** で足りる |
| `C ≥ 2` の由来 | `(1/2)(β/α)² ≥ β/α` の**ちょうどの条件** |

### ★「正面から要ると思ったものが要らなかった」が 3 件増えた(累計 11 件)

9. **単調性に導関数は要らなかった** —— 性質 (e) から直接出る
10. **`f₀ ≤ 1/4` に 4 場合分けは要らなかった** —— AM–GM 1 回
11. **性質 (d) に `n` の偶奇の場合分けは要らなかった** ——
    `f₀ = max(0, −f(r))` が両方の場合を 1 つの式で表す

### ★自己訂正 1 件(本日 3 度目)

`IsBelyiPoly.comp`(Belyi 写像は合成で閉じる)を
「`Lemma 2.2` の帰納法の構造的な核」と書いたのは**誤りだった**。
帰納法が合成するのは `g ∘ (f + f₀)` で、`f + f₀` は `0, 1` を `f₀` へ写すので
**Belyi 型ではない**。合成が閉じるのは `f + f₀` の臨界値が `f(S) ⊆ S′` に入るからで、
相対版 `comp_crit_of_rel` がそれを取る。

★本日の自己訂正 3 件はいずれも**その日のうちに**取り消せた。
**記録の粒度が細かいほど誤りが早く出る。**

## §9-9 2026-08-17 —— §1 の律速を測り直した:**複素解析空間は要らない**

### ★★まず、`Proposition 1.4` の statement は**偽**だった

`Skeleton` は `∀ D : HeightTheoryData, …` と量化しているが、
**`HeightTheoryData` は公理を 1 つも持たないデータ**である。
★反例を `Check/GenEll/HeightAxiomGap.lean` に構成した(機械検証済み)——
`Point ≝ ℕ`、`ABundle ≝ ℝ`、`tensor ≝ (+)`、`ht L x ≝ L²` で
**(i) 加法性と (iv) Northcott が同時に破れる**。

★★**あの `sorry` は「まだ証明していない」ではなく「証明できない」だった。**
★`check.mjs` 冒頭 B5(条件を posit して `sorry` を消す穴)の**裏側**である——
**条件を posit しないまま全称量化すると statement が偽になる。**

★★閉じる道は **`HeightTheoryData` を posit ではなく構成に置き換える**しかない。
公理を足すのは B5 そのものである。

### ★★★`Definition 1.1` を層に分けた —— 複素解析空間は要らない

原典 `ht` の定義:
> ★★この 1 行が層の境界を決めている——一般の `X` 上の算術直線束が要るのは、
> `x_F^* M̄ ∈ APic(Spec 𝒪_F)` を作るためだけである。

★`x_F^* M̄` の中身は:
- **有限側**: `x_F^* L` = 可逆 `𝒪_F` 加群 = 分数イデアル類。**スキーム論だけ。解析は要らない**
- **アルキメデス側**: 各 `v : InfinitePlace F` について、**1 点** `x_v ∈ X(ℂ)` の
  ファイバー上のノルム。★★**連続性も位相も要らない**

★★★**したがって `ht` の定義そのものには、位相も解析的構造も要らない。**

| 項目 | 要るもの |
|---|---|
| `Definition 1.1` のデータ / `ht` / `Definition 1.2, (i)` | ★**可逆層 + 点ごとのノルム**のみ |
| `Proposition 1.4, (i)`(加法性) | 同上(テンソルのノルムは積) |
| `Proposition 1.4, (ii)(iii)` / `Proposition 1.6` のアルキメデス側 | ★**`X(ℂ)` の位相 + コンパクト性 + ノルムの連続性** |
| `Proposition 1.4, (iv)`(Northcott) | ample + 射影埋め込み → ★`NorthcottProj.lean` は実装済み |
| `Example 1.3, (ii)`(compactly bounded) | `X(ℂ_v)` の位相(非アルキメデス側も) |

★★**「複素解析空間」——正則構造・連接層・GAGA——は §1 のどこにも現れない。**
現れるのは**位相とコンパクト性と連続性**だけである。

### ★mathlib 実測(2026-08-17)

| 要るもの | mathlib |
|---|---|
| スキーム上の**可逆層 / 直線束** | ★**無い**(`AlgebraicGeometry/Modules/` は `Presheaf`/`Sheaf`/`Tilde` の 3 本) |
| `Projectivization` の**位相** | ★**無い**(instance は代数・群作用のみ) |
| `ℙⁿ(ℂ)` の**コンパクト性** | ★**無い** |
| 複素解析空間 / GAGA | 無い(★**しかし要らない**) |

### ★★★規模の見直し

昨日までの記録は `Definition 1.1` を「`X^arc`(複素解析空間)待ち = 研究規模」としていた。
★**訂正する。** 要るのは

1. **スキーム上の可逆層**(引き戻しつき)——mathlib 規模の仕事
2. **`ℙⁿ(ℂ)` の位相とコンパクト性**——商位相 + 単位球面の像。数百行
3. 閉部分スキームの `ℂ` 値点が閉(→ コンパクト)——2 の系

★★**研究規模ではなく、数週間〜数か月規模である。**
★これは本日 3 度出た「名前で測って『無い』と結論した」誤りと同じ型の訂正である
(`Cartier` / `scheme 論的像` に続く 4 例目)。

## §9-10 2026-08-17 —— 可逆層は正面から作らない:**Arakelov 因子の提示**

### ★★mathlib 実測: `SheafOfModules` にテンソル積も引き戻しも無い

| 要るもの | mathlib |
|---|---|
| `X.Modules = SheafOfModules X.ringCatSheaf`(圏・極限・余極限・自由加群) | ★**ある** |
| **テンソル積**(モノイダル構造) | ★★**無い**(`tensorObj` / `MonoidalCategory` で 0 件) |
| **引き戻し `f^*`** | ★★**無い**(`ChangeOfRings.lean` は `restrictScalars` のみで逆向き) |

★★**したがって可逆層を正面から作ると、先に「加群層のモノイダル構造」
(前層のテンソル積の層化)を作ることになる——mathlib 規模の仕事である。**

### ★★★迂回路 —— Cartier 因子で提示する

正規整スキーム `X` について **`Pic(X) ≅ CaCl(X)`**(Cartier 因子を主因子で割ったもの)。
★**この提示を取れば層のテンソル積は要らない**:

| 演算 | 層で書くと | 因子で書くと |
|---|---|---|
| `L ⊗ M` | ★モノイダル構造(無い) | ★**因子の和**(自明) |
| `x_F^* L` | ★層の引き戻し(無い) | ★**因子の引き戻し**——`CartierPullback.lean` に**実装済み** |
| `Pic(X)` | 同型類 | ★**因子 / 主因子** |

★★★**これは `ADiv(F)` と同じ形である**——有限側が因子、アルキメデス側が実数。
`APic(X) ≝ (Cartier 因子) × (X(ℂ) 上の Green 関数) / (主因子)`。

### ★既に手元にある部品

| 部品 | 場所 |
|---|---|
| `IsEffectiveCartier` / `pullbackIdeal` / `conductorADiv` | `Found/GenEll/CartierPullback.lean` |
| `ADiv F` / `deg` / `degNormalized` / `APicOF F` | `ArithDiv.lean` / `ProductFormula.lean` |
| `degNormalized` の底変換不変性 | `BaseChange.lean` |
| **`ℙ(V)` の位相とコンパクト性** | `ProjTopology.lean` |
| 射影空間の Northcott / 高さの拡大公式 / 絶対高さ | `NorthcottProj.lean` / `HeightExtension.lean` |

★★**`Definition 1.1` を構成するのに新たに要るのは:**

1. Cartier 因子の**群**(有効なものの差)—— `IdealSheafData` の上で
2. **主因子**(有理関数から)
3. **Green 関数**(= `X(ℂ)` 上の計量。`-log‖s_D‖` の形)
4. `x_F` による引き戻し(1〜3 を通す)—— 有限側は**実装済み**

### ★これが本日 5 例目の「正面は要らなかった」になる可能性

`Cartier` / scheme 論的像 / NCBelyi の逐語 / 複素解析空間 に続く。
★★**まだ確認していない**——上の 1〜4 を実際に取るまでは仮説である。

---

## §9-11 2026-08-17 夕 —— **迂回路は通った**:高さを構成し、`Prop 1.4` の (ii)(iii) を無条件で取った

§9-10 は「`Definition 1.1` を Cartier 因子で提示する」という迂回路を提案し、
「まだ確認していない——1〜4 を実際に取るまでは仮説である」と書いた。
★★★**同日中に 1〜4 のうち 3 つを取り、高さそのものを構成した。**

### ★★★何を取ったか

| §9-10 の項目 | 状態 | 場所 |
|---|---|---|
| 1. Cartier 因子のモノイド構造 | ★★★**取得** | `IdealStalk.lean`(`isEffectiveCartierStalk_mul`) |
| 2. 主因子 | ★`Spec 𝓞_F` 側は取得済(`APrc`) | `ArithDiv.lean` / `ProductFormula.lean` |
| 3. Green 関数 | ★★★**取得** | `ArchPoint.lean`(`GreenFn` / `archADiv`) |
| 4. `x_F` による引き戻し | ★★★**取得**(有限側+無限側) | `ArchPoint.lean`(`pullbackADiv`) |
| **高さそのもの** | ★★★**構成した** | `HeightConstruction.lean`(`htArith`) |

### ★★★`Definition 1.1` が要求する `X^arc` は、引き戻しには要らなかった

原文は「コンパクト正規複素解析空間 `X^arc` 上のエルミート計量」を要求する。
★★引き戻し `x_F^* L̄` に実際に要るのは **ℂ-点だけ**である——
`𝓞_F → F ↪_v ℂ` の `Spec` を取って `x_F` に合成すればよい。

★★★**これが「正面から要ると思ったものが要らなかった」の 6 例目**である
(`Cartier` / scheme 論的像 / NCBelyi の逐語 / 複素解析空間 / 茎による定義 に続く)。

### ★★★`Proposition 1.4` の 4 条の現状

| 条 | 構成側の状態 | 仮説 |
|---|---|---|
| (i) 加法性 | 取得 | ★`PullbackMul`(= `comap` の積保存)1 本 |
| (ii) 非負性 | ★★**無条件で取得** | 無し(Green 関数の非負性のみ) |
| (iii) 類にしか依らない | ★★**無条件で取得** | 無し |
| (iv) Northcott | 未着手 | —— |

★★**(ii) は原文より強い形で出た。** 原文の結論は `≳ 0`(定数差を許す)だが、
構成では **`0 ≤ ht`(定数を要しない真の不等式)**が出る。
`Definition 1.2, (ii)` の**印字の向きの食い違い**は、ここでは争点にならない。

★★**(iii) も原文より強い。** 原文は BD-class の一致だが、
構成では **`ht` そのものが一致する**(定数差すら出ない)。
機構は `deg` が主算術因子の上で消えること(積公式)であり、
`degNormalizedAPic : APic(Spec 𝓞_F) →+ ℝ` は既に構成済みだった。

### ★★★原文の仮定の所在が変わった —— (ii) の場合

原文 (ii) は「`L_ℚ` のある正冪が大域切断で生成される」を仮定する。
★★因子で表すと、有限素点側は**イデアル層**(`O_X` の部分層)であり、
引き戻しても `𝓞_F` の**イデアル**であって分数イデアルにならない。
★★★**つまり原文の仮定は「有効因子で代表できる」ことに吸収されており、
残る仮定はグリーン関数の非負性だけである。**

### ★★★§1 の `Skeleton` は**系統的に**偽である(3 件を機械検証)

`Check/GenEll/HeightAxiomGap.lean`(既存)に加え、
`Check/GenEll/RemarkAxiomGap.lean` を新設して
**`Remark 1.4.1` と `Remark 1.5.1` も偽である**ことを機械検証した。

★反例には新しい道具が要った——既存の `badData` は `ht L x = L²` で
**`x` に依らない**ので、定数差を許す `BDeq` は成り立ってしまう。
★★**点に線形に依る高さ**(`linData`: `ht L x ≝ L · x`)が必要だった。

★★★**これは `Proposition 1.4` 1 件の問題ではなく、
`HeightTheoryData` を posit したことの帰結である。**
閉じる道は 1 つ——posit を構成に置き換えることであり、本日それを実行した。

### ★★`Definition 1.5, (ii)` の残りが 1 本に絞れた

`RadicalCartier.lean` で茎と根基の交換・点ごとの主張・scheme への大域化・
被約性の冪等性をすべて取得した。
★★★**残るのは Auslander–Buchsbaum(正則局所環 ⇒ UFD)ただ 1 本**である。

### ★★★器具の記録 —— 同じ数学が「書き方」で通ったり通らなかったりする(本日 3 例)

| 件 | 通らない書き方 | 通る書き方 |
|---|---|---|
| `germ_res` | `rw` / `simp only` | ★**`erw`**(`X.presheaf` が 2 通りに解釈される) |
| `Ideal.comap_symm` | `rw`(`↑e.symm` と `e.symm` で強制の経路が違う) | ★**`have := fun J => …` で定義的等しさで渡す** |
| `IsLocalization.map_radical` | `⟨x, hxU⟩` を組む(**`whnf` が 100 万ヒートビートでも終わらない**) | ★★**`y : ↥U` を変数として持つ補題を経由** |

★★★**3 件とも数学は 1 ミリも変わっていない。**
★これは「実装の摩擦」ではなく**測定すべき現象**である——
mathlib の型階層と強制の経路が、証明の書き方を決めている。

### ★器具が私の誤りを 2 件捕まえた

1. **逐語の捏造**: `> (i) We have: ht` と書いたが、原文 p.6 の該当行は
   `(i) We have` のみで、**コロンも `ht` も無い**(数式は別行に組まれている)。
   `check.mjs` の引用照合が NG 2 件で落とした。
2. **`sectionId` の捏造**: `genell-def-1-5-ii` と書いたが、
   構造化 HTML にあるのは `genell-def-1-5` だけである。
   ★**同じ内容の `RadicalPrincipal.lean` が既に正しい id を使っていたのに、
   確認せずに書いた。** G1 が NG 2 件で落とした。

### ★★残る穴は「mathlib に無い 1 本の定理」に収束した

| 項目 | 残り | 種別 |
|---|---|---|
| `Prop 1.4, (i)` | `Scheme.IdealSheafData.comap_mul` | ★mathlib 級の貢献 1 本 |
| `Def 1.5, (ii)` | Auslander–Buchsbaum | ★★mathlib に着手の形跡なし |

★★★**どちらも器具の問題ではなく、数学の在庫の問題である。**
★`blocked-leaves.json` に両方を記録した(後者は既存 entry を更新、前者は新規)。

### ★進捗の数字は動いていない(5/24)

★**これは正しい状態である。** 1 件に数えるには命題**全体**が要り、
`Proposition 1.4` は (i)〜(iv) の 4 条すべてが揃って初めて 1 件になる。
★★条つき `.src` を数に入れれば数字は上がるが、それは B5 そのものである。

### ★§9-11 追記 —— 迂回路の部品を精密に測った(2026-08-17 夜)

`comap_mul` を迂回する道(台の有限性 + 茎での積保存)の部品を、
実際に `lake env lean` で当てて測った。

| 部品 | 状態 |
|---|---|
| 有限素点 → `Spec 𝓞_F` の点 | ★**そのまま通る**(`⟨v.asIdeal, v.isPrime⟩`) |
| 茎 ≅ 局所化 | ★`StructureSheaf.stalkIso R p` は **`AlgEquiv`** で在る |
| `Algebra 𝓞_F (茎)` の instance | ★★**無い**——`stalkIso` を明示的に通す必要がある |
| 局所化が DVR | ★`IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain` が在る |
| 　　同 instance | ★★**自動では出ない**(`P ≠ ⊥` と `IsLocalization.AtPrime` を明示する) |
| 台が有限 | ★★`support_comap`(mathlib)+ `stalk_eq_top_iff`(本日取得)で道は通っている |
| 茎での積保存 | ★★★**取得済**(`stalkPullback_mul`) |

★★**数学の穴は無い。**すべての部品が mathlib か本プロジェクトに在る。
★残るのは instance を明示的に通す配線であり、
**`Localization.AtPrime v.asIdeal` の側で作業して `stalkIso` で最後に移す**のが筋である
(茎の側で instance を合わせようとすると `RadicalCartier.lean` と同じ摩擦に入る)。

★★★**次のセッションはここから**——半端に作らず、測定だけを残す。

---

## §9-12 2026-08-17 夜 —— **`X^arc` の判定を撤回した**:複素解析空間は要らない

### ★★★2026-08-16 の判定は誤りだった

`blocked-leaves.json` の `[GenEll] Definition 1.1` は
「複素解析空間の理論が要る。mathlib に着手の形跡が無い」と書いていた。
★★★**撤回する。** `Definition 1.1` が実際に要求するのは
`X^arc` の**位相とコンパクト性**だけであり、
正則構造・連接層・GAGA はどこにも現れない。

### ★★本日 3 段階で解いた

| 段 | ファイル | 使ったもの |
|---|---|---|
| 1. `ℙ(V)` がコンパクト | `ProjTopology.lean` | 商位相 + 単位球面の連続像。★Hausdorff 性は不要 |
| 2. 射影集合が閉 | `ProjClosed.lean` | ★**多項式の連続性だけ**。GAGA 不要 |
| 3. `X^arc` に位相が載る | `ArcModel.lean` | 射影モデルをデータとして受け、コンパクト性を**証明** |

★さらに `ArchPoint.lean` で **ℂ-点(`Spec ℂ ⟶ X`)を数論的な埋め込みから構成**し、
Green 関数の引き戻しと `htArith`(高さそのもの)まで作った。

### ★★★`ArcModel` は posit ではない

`Interface` の `structure` は「まだ作っていないもの」を受ける。
★★`ArcModel` が受けるのは**埋め込みそのもの**——
「与えられた `X` について実際に構成できるもの」であり、**未証明の事実ではない**。
コンパクト性は `ArcModel.compactSpace` で**証明した**。

★★★**この区別が B5 を避ける線である。**
posit は「事実を仮定する」、データは「対象を受け取る」。

### ★★`Definition 1.1` に残るのは 2 本

| 残り | 種別 |
|---|---|
| 射影埋め込みの構成(`ArcModel` を実際に作る) | ★代数幾何の配線 |
| `Scheme.IdealSheafData.comap_mul` | ★mathlib 級の貢献 1 本 |

★★**どちらも複素解析ではない。**

### ★`Example 1.3, (ii)` も posit を外せる

原文の「`x` が定める `X^arc` の `[F:ℚ]` 個の点」は `archPoint` そのものだった。
★`CompactBound.lean` で `BoundedByArch` を**定義**し、
「囲われた点ではアルキメデス側の寄与が一様に抑えられる」ことを**証明**した。
★★★定数が `F` に依らないのは**正規化のおかげ**である
(`#{無限素点} ≤ [F:ℚ]`)。

### ★本日の `Found/GenEll/` 新規ファイル(14 件、すべて sorry 0・標準 3 公理のみ)

```
IdealStalk  StalkPullback  DegMul  ArchPoint  HeightConstruction
HeightNonneg  HeightClass  RadicalCartier  StalkSupport
ConductorHeight  ArchBound  ProjClosed  CompactBound  ArcModel
```

### ★進捗の数字は 5/24 のまま

★**これは正しい状態である。** 1 件に数えるには命題**全体**が
原文どおりの仮定で揃う必要がある。
★★条つき `.src` を数に入れれば数字は上がるが、それは B5 そのものである。

---

## §9-13 2026-08-17 深夜 —— **器具の規則が 1 つ出た**

### ★★★規則: 圏の言葉で終わらせ、元には最後に 1 度だけ降りる

本日 5 度、**同じ数学が書き方で通ったり通らなかったり**した。
そのうち最後の 1 件で、**対処法が言語化できた**。

| 件 | 通らない書き方 | 通る書き方 |
|---|---|---|
| `germ_res` | `rw` / `simp only` | `erw` |
| `Ideal.comap_symm` | `rw`(強制の経路が違う) | `have := fun J => …` で定義的等しさで渡す |
| `IsLocalization.map_radical` | `⟨x, hxU⟩` を組む(`whnf` 発散) | `y : ↥U` を**変数として持つ**補題を経由 |
| `isClosed_coinduced` | `rw`(`inferInstanceAs` が構文的に合わない) | `convert … using 1` |
| `comap` の易しい包含 | 途中で `simp only` を当てる | ★★★**射のレベルで分解を完成させ、最後に 1 度だけ降りる** |

### ★★摩擦の正体

`CommRingCat.Hom.hom` と `ConcreteCategory.hom` の**正規形が揺れ戻る**。
`simp only` を当てるたびにどちらかへ寄り、次の `rw` が外れる。
`CommRingCat.comp_apply` / `ConcreteCategory.comp_apply` / `RingHom.comp_apply` /
`CommRingCat.hom_comp` のどれを入れても別の形になる。

★★★**同じ証明が import の重さで通ったり通らなかったりした**——
軽い import の probe では通り、本ファイルに置くと通らない。
★`simp only` は与えた補題しか使わないはずだが、
**環境にある instance が正規化の経路を変える**。

### ★★★対処

1. **証明を圏の言葉(`≫`、`Category.assoc`、`Scheme.Hom.comp_app`)で完成させる**
2. `eqToHom` の証明項が違って `rw` が照合に失敗する所は
   **`exact Category.assoc _ _ _`** で当てる(証明無関係で通る)
3. **最後に 1 度だけ** `CommRingCat.hom_comp` + `RingHom.comp_apply` で元に降りる
4. ★★**途中で `simp` を 1 度も使わない**

★これで `map_appLE_le_ideal_comap`(10 回失敗していた)が通った。

### ★この規則が効く範囲

`Found/FrdI/` を含め、**圏論を経由するすべての証明**に効く。
★★「実装の摩擦」として流さず、**規則として書いた**ことが本日の成果である。

---

## §9-14 2026-08-17 深夜 —— **原文が要求する「複素共役との両立」の理由が出た**

### ★★★`Definition 1.1` はなぜ `ι_X` との両立を要求するのか

原文 `Definition 1.1` は、エルミート計量が
**複素共役 `ι_X` と両立する**ことを要求する。★原文は理由を書いていない。

★★★**理由は「高さの底変換不変性」である。**

### ★★機構

高さが `X(ℚ̄) → ℝ` として well-defined であるには、定義体 `F ⊆ K` の
取り替えで値が変わらないこと(底変換不変性)が要る。

- **有限素点側**: `x_K^* D = (x_F^* D)·𝓞_K`
  ★本日 `PullbackBase.lean` の `pullbackIdeal_comp` で取得
- **アルキメデス側**: `w : InfinitePlace K` を `F` に制限すると `v` になるが、
  ★★**`v.embedding` は `w.embedding|_F` **または その複素共役**に等しい**
  (mathlib `NumberField.InfinitePlace.embedding_mk_eq`:
  `embedding (mk φ) = φ ∨ embedding (mk φ) = conjugate φ`)

★★★**したがって Green 関数が複素共役で不変でなければ、
`archPoint` の値が `F` の取り方に依ってしまう。**

### ★★これで「原文が書いていない理由」の 2 例目

| 原文の要求 | 書かれていない理由 | 見つけた場所 |
|---|---|---|
| 次数を `[F:ℚ]` で**正規化**する | ★`≲` の定数を `F` に依らなくするため | §9-11(`ArchBound.lean`) |
| 計量が `ι_X` と**両立**する | ★★高さの**底変換不変性**のため | ★本節 |

★★★**どちらも「型が要求したので探したら、原文の中に既にあった」**という形である。
Gap(原文の飛躍)ではない——**原文は正しく、理由を省いていただけ**である。

### ★★`ArithCartier` に条件を足す必要がある

`Found/GenEll/ArchPoint.lean` の `ArithCartier` は
Green 関数に条件を付けていない。
★★★**底変換不変性を取るには `green` に複素共役不変性を課す必要がある。**
それは `Definition 1.1` の忠実な写しでもある。

★次のセッションはここから: `ArithCartier` に `green_conj` を足し、
アルキメデス側の底変換を取り、`degNormalized_baseChange` と繋ぐ。
★それが入れば `ht` が `X(ℚ̄)` の上で well-defined になり、`Definition 1.2` が閉じる。

---

## §9-15 2026-08-17 深夜 —— `Definition 1.2` に残る鎖を測った

### ★★底変換の両側は取れた

| 側 | 場所 | 状態 |
|---|---|---|
| 有限素点側 `x_K^* D = (x_F^* D)·𝓞_K` | `PullbackBase.lean` | ★★★取得 |
| アルキメデス側(ℂ-点の値) | `ArchBaseChange.lean` | ★★★取得 |

★アルキメデス側は共役の場合分けを `IsConjInvariant`(原文の `ι_X` 両立)が吸収する。

### ★★★残る鎖(次のセッションはここから)

1. **有限素点側を `ADiv` に繋ぐ**:
   `idealADiv K (J·𝓞_K) = baseChangeFin F K (idealADiv F J)`
   ★要るのは **`ord_w(J·𝓞_K) = e(w/v)·ord_v(J)`**。
   ★★mathlib には `ramificationIdx`(素イデアルについて)はあるが、
   **一般のイデアルの拡大の重複度公式は無い**(2026-08-17 実測)。
   分解の乗法性から組む必要がある。

2. **アルキメデス側を `ADiv` に繋ぐ**:
   `archADiv K g xK = baseChangeArc F K (archADiv F g xF)`
   ★要るのは無限素点の相対次数和(`InfinitePlaceRel.lean` に部品あり)。

3. **`degNormalized_baseChange` を当てる**(既に取得済)。

4. **`X(ℚ̄)` の型を作る**: 数体の圏についての colimit。
   ★★これは設計の仕事であり、証明の仕事ではない。

### ★★これが取れると何が起きるか

`ht : X(ℚ̄) → ℝ` が well-defined になり、
`Definition 1.2, (i)` が閉じる。(ii) は既に取得済(`BDClass.lean`)なので
★★★**`Definition 1.2` が完成し、指標が 5/24 → 6/24 に動く。**

★本日の指標が動かなかったのは、この 4 段が残っているためである。
★★**数学の穴は 1 だけ**(一般のイデアルの拡大の重複度公式)であり、
2 は部品あり、3 は取得済、4 は設計である。

### ★★★§9-15 の即時訂正 —— 「数学の穴」は無かった

★直前に「一般のイデアルの拡大の重複度公式は mathlib に無い(実測)」と書いた。
★★★**誤りである。取り消す。**

    `NumberField/RamificationInertia/Ramification.lean:310`
    theorem emultiplicity_map_eq_ramificationIdx_mul
      [IsDedekindDomain R] [FaithfulSMul R S]
      {v : Ideal R} {w : Ideal S} {I : Ideal R} (h : I ≠ ⊥)
      (hv : Irreducible v) (hw : Irreducible w) (hw_bot : w ≠ ⊥) [w.LiesOver v] :
      emultiplicity w (I.map (algebraMap R S))
        = v.ramificationIdx w * emultiplicity v I

★★**`I` は一般のイデアルである。** 私は `ramificationIdx` を素イデアルの語で検索し、
`emultiplicity` の語で検索しなかった。

★★★**測り方の誤り**である——「名前で測って『無い』と結論した」の再発。
本日の午前に `Cartier` で同じ誤りをして訂正したばかりである。

### ★★これで `Definition 1.2` に数学の穴は無い

| 段 | 状態 |
|---|---|
| 1. 有限素点側を `ADiv` に繋ぐ | ★**部品は mathlib にあった**(`emultiplicity_map_…`)。`Associates.count` との橋が要る |
| 2. アルキメデス側を `ADiv` に繋ぐ | 部品は `InfinitePlaceRel.lean` |
| 3. `degNormalized_baseChange` | 取得済 |
| 4. `X(ℚ̄)` の型 | 設計 |

★★★**4 段すべてが配線と設計であり、未知の数学は無い。**

### ★段 1 の橋も測った(2026-08-17 深夜)

`idealADiv` は `Associates.count _ (Associates.mk J).factors` を使い、
mathlib の重複度公式は `emultiplicity` で述べられている。
★**両者を直接繋ぐ補題は無い**(`UniqueFactorizationDomain/` を実測)。

★★ただし**どちらも整除で特徴づけられる**ので橋は架かる:

- `Associates.prime_pow_dvd_iff_le` —— `p^n ∣ a ↔ n ≤ count p a.factors`
- `emultiplicity` の定義 —— `p^n ∣ b ∧ ¬p^(n+1) ∣ b`

★★★**したがって段 1 も「未知の数学」ではなく配線である。**
次のセッションはこの橋から始めるのが最短である。

---

## §9-16 2026-08-17 深夜 —— `Definition 1.2` の証明は終わり、**設計の問題が 1 つ残った**

### ★★★仮定はすべて定理になった

`htArith_baseChange`(高さの底変換不変性)が受けていた仮定は、
本日すべて定理になった:

| 仮定 | 定理にした場所 |
|---|---|
| `hfin`(有限素点側) | `PullbackNatural.lean`(`ΓSpecIso` の自然性) |
| `harch` の埋め込み両立 | `ArchCompat.lean`(タワー + 共役の場合分け) |
| `hlies` | `FinitePlaceRel.lean`(instance) |

★★★**未証明の数学は 1 つも残っていない。**

### ★★★残った設計の問題 —— `ArithCartier` は**モノイドであって群でない**

`htArith_baseChange` はもう 1 つ `hJ : x_F^* D ≠ 0` を要求する。
★これは「`x` が `D` を通らない」という条件であり、**一般には成り立たない**。

★★原因は `ArithCartier` の `divisor` が**有効 Cartier 因子(イデアル層)**だからである。
原文の算術直線束は `APic(X)` の元、すなわち**因子の差**である。

★★★**したがって `X(ℚ̄) → ℝ` を無条件に定めるには、
`ArithCartier` を有効因子の**差**に拡張する必要がある。**

- 有効因子 2 つの組 `(D⁺, D⁻)` として持つ
- `x_F^*` は成分ごと、`deg` は差
- ★`x` が `D⁺` にも `D⁻` にも通らない、という条件は依然要るが、
  **一般の直線束は「通らない」代表を選べる**(移動補題)

### ★★これで `Definition 1.2` に残るもの

| 項目 | 種別 |
|---|---|
| `ArithCartier` を差に拡張 | ★設計(証明は既存のものが成分ごとに効く) |
| `X(ℚ̄)` の型(数体の colimit) | ★設計 |
| 移動補題(代表の選択) | ★★数学。ただし `Definition 1.2` の**名前づけ**には要らない |

★★★**証明の仕事は終わった。次のセッションは設計から始める。**

### ★★★§9-16 の追測 —— 差にしても `hJ` は消えない。**それが正しい**

「`ArithCartier` を差に拡張すれば `hJ` が消える」と書いたが、★**消えない**。
差の形でも各成分に `x_F^* D^± ≠ 0` が要る。

★★★**しかしそれは欠陥ではなく、因子表示の正しい境界である。**

原文の `ht_L̄(x) = deg_F(x_F^* L̄)` は**可逆層**の引き戻しを使う——
可逆層はいつでも可逆層に引き戻るので、条件が要らない。
★★因子で表すと「`x` が `D` を通らない」が要るのは、
**表示の側の制約**であって高さの側の制約ではない。

### ★★これが意味すること

| 対象 | 因子表示での状態 |
|---|---|
| `ht` を `X(ℚ̄) \ D` の上で定める | ★★★**本日すべて取れた**(仮定は全部定理) |
| `ht` を `X(ℚ̄)` **全体**で定める | ★移動補題(`D` を動かして `x` を避ける)が要る |

★★★**原文が `Proposition 1.6` を `U(ℚ̄) = X(ℚ̄) \ D` の上で述べているのは、
まさにこの境界に沿っている。**
本日構成した `htArith` は、原文が実際に使う場所では**そのまま使える**。

### ★次のセッションの選択肢は 2 つ

1. **移動補題**を作る(`D` を線形同値で動かす)—— 数学の仕事
2. **可逆層**で表す —— mathlib に `SheafOfModules` のテンソル積が無いので迂回が要る

★★どちらも `Definition 1.2` の**名前づけ**(`ht_M̄ : X(ℚ̄) → ℝ` と呼ぶこと)には
必要だが、★**`Proposition 1.6` や `Theorem 2.1` が使う範囲では要らない**。

---

## §9-17 2026-08-17 深夜 —— **`Definition 1.2, (i)` が `U_X(ℚ̄)` の上で完成した**

### ★★★到達点

    `htU : U_X(ℚ̄) → ℝ`   (`Found/GenEll/UPoint.lean`)

原文の高さ関数そのものである。★仮定は原文自身のもの 2 つだけ:

- `IsConjInvariant D.green` —— 原文の「計量は `ι_X` と両立する」(`Definition 1.1`)
- 各点の `x_F^* D ≠ 0` —— 型 `AlgPointOff` に組み込まれている

### ★★★積み上げた鎖(すべて本日、すべて sorry 0)

| 段 | ファイル |
|---|---|
| ℂ-点と Green 関数 | `ArchPoint.lean` |
| 高さの構成 | `HeightConstruction.lean` |
| `comap` は積を保つ | `ComapLocal` → `ComapAffine` → `ComapMul` |
| `Prop 1.4, (i)` が無条件に | `HeightAdditive.lean` |
| 引き戻しの底変換(有限側) | `PullbackBase` → `PullbackNatural` |
| 複素共役との両立 | `ArchConj` → `ArchBaseChange` → `ArchCompat` |
| `ADiv` への翻訳 | `IdealADivBase` / `ArchADivBase` |
| 高さの底変換不変性 | `HeightBaseChange` → `HeightInvariant` |
| 型と商 | `AlgPoint` → `UPoint` |

### ★★なぜ指標が動かないのか(正直に)

原文の `Definition 1.2, (i)` は **`X(ℚ̄)` 全体**で高さを定める。
★因子表示では「`x` が `D` を通らない」が要るので、
取れたのは **`U_X(ℚ̄) = X(ℚ̄) \ D`** の上である。

★★全域化には 2 つの道があり、**どちらも別規模**である:

1. **可逆層で表す** —— 可逆層はいつでも引き戻せるので条件が要らない。
   ★mathlib に `SheafOfModules` のテンソル積が無いので迂回が要る
2. **移動補題** —— `D` を線形同値で動かして `x` を避ける。
   ★主因子(`X` 上の有理関数)の理論が要る

★★★**`.src` は条つきに留めた。** 条なしにすれば 5/24 → 6/24 に動くが、
**原文が定める範囲を狭めて数えるのは B5 である。**

### ★原文が実際に使う範囲では完成している

`Proposition 1.6` は `U_X(ℚ̄)` の上で述べられている。
`Theorem 2.1` は compactly bounded subset の上で述べられている。
★★**どちらも本日の構成がそのまま使える。**
