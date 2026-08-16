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
