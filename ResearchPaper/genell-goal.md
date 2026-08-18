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

---

## §9-18 2026-08-17 深夜 —— `Proposition 1.6` が組み上がり、残る 1 点を測った

### ★★★`Proposition 1.6` は組み上がった

    `log-cond_D(x) ≤ ht_D(x) + C`   (`C` は `F` にも `x` にも依らない)

`Found/GenEll/Prop16.lean`(sorry 0、標準 3 公理のみ)。
★**射影モデルと Green 関数の連続性さえ与えれば、定数は自動で出る。**

★★`L = O_X(D)` という原文の対応は、因子表示では**定義そのもの**になった。

### ★★★残るのは `ArcModel` の構成 1 点

`ArcModel`(`X(ℂ)` を `ℙ(V)` の閉錐集合として実現するデータ)を
**与えられたものとして受けている**。
★原文は「`X` は ℤ-固有」から射影埋め込みを得る
(`Remark 1.4.1` が「ℤ-projective model that arises from a projective embedding」と明記)。

★★構成に必要な段:

| 段 | 状態 |
|---|---|
| `X(ℂ) → ℙⁿ(ℂ)`(閉埋め込みとの合成) | ★自明 |
| 単射性(閉埋め込みは mono) | ★mathlib にある |
| **`ℙⁿ(ℂ) ≅ ℙ ℂ (Fin (n+1) → ℂ)`** | ★★★**mathlib に無い** |
| 像が閉(方程式で切れる) | ★`ProjClosed.lean` で取得済 |

### ★★★測定: `Proj` の「点の関手」は mathlib に無い

`Mathlib/AlgebraicGeometry/ProjectiveSpectrum/` は
`Basic` / `Functor` / `Proper` / `Scheme` / `StructureSheaf` / `Topology` の 6 ファイル。
★`Functor.lean` にあるのは **`Proj` の間の射の関手性**であり、
**「`Proj` への射 = 直線」という点の関手の記述は無い**(2026-08-17 実測)。

★★これは mathlib 級の貢献 1 本である。

### ★★これで §1 の残りが出そろった

| 項目 | 残り | 種別 |
|---|---|---|
| `Prop 1.6` | `ℙⁿ` の点の関手 | ★★mathlib 級 1 本 |
| `Def 1.2, (i)` | 可逆層 or 移動補題 | ★★別規模 |
| `Def 1.5, (ii)` | Auslander–Buchsbaum | ★★mathlib に着手なし |
| `Prop 1.4, (iv)` | Northcott(ample との比較) | ★★別規模 |

★★★**どれも「原文の数学」ではなく「mathlib の在庫」の問題である。**

### ★★★`ℙⁿ` の点の関手 —— 構成の 5 段を測った(2026-08-17 深夜)

`Proposition 1.6` と `Proposition 1.4, (ii)` の両方が、
`ArcModel`(`X(ℂ)` を `ℙ(V)` の閉錐集合として実現するデータ)の**構成**に掛かっている。

★mathlib 実測: `Proj` の**アフィン被覆はある**が、**点の関手は無い**。

| 部品 | 状態 |
|---|---|
| `Proj.awayι : Spec (Away 𝒜 f)₀ ⟶ Proj 𝒜` | ★在る(標準チャート) |
| `Proj.basicOpen` / `opensRange_awayι` | ★在る |
| `Proj.affineOpenCover` | ★在る |
| **「`Proj` への射 = 直線」という点の関手** | ★★★**無い** |

### ★★構成の 5 段

1. `p : Spec ℂ ⟶ Proj 𝒜` の像が `Proj.basicOpen x_i` に入る `i` を取る
   (★無関係イデアルを含まないから存在する)
2. 開埋め込みなので `p` は `awayι` を経由して分解する
3. 環準同型 `(Away 𝒜 x_i)₀ → ℂ` を得る——これが座標 `x_j/x_i ↦ c_j` である
4. 射影点 `[c_0 : … : 1 : … : c_n] ∈ ℙ ℂ (Fin (n+1) → ℂ)` を作る
5. `i` の取り方に依らないこと・単射性を示す

★★★**これは mathlib 級の貢献 1 本である**(`HomogeneousLocalization` の摩擦つきで
200 行規模)。★取れれば **2 件同時に**動く(`Prop 1.6` と `Prop 1.4, (ii)`)。

### ★なぜ `[CompactSpace (complexPoints X)]` を仮定にしないのか

そうすれば両命題は「条なし」に見える。★★**しかしそれは B5 である。**
原文の `Definition 1.1` は `X^arc` が**コンパクトであること**を
`X` が ℤ-固有であることの帰結として述べており、**仮定ではない**。

★★★`ArcModel` は**データ**(埋め込み)を受けてコンパクト性を**証明する**——
そこが posit との違いである。

---

## §9-19 ★★★残作業を条ごとに測り直した——記録していた見積りが 2 つ誤っていた

**2026-08-17 深夜。**カウンタが動かない理由を「怠慢ではなく 2 値規則」と書いてきたが、
★**その先の「あと何が要るか」を条の水準で測っていなかった**。測り直した。

### ★★★分子の規則(再確認)

`tools/genell-progress.mjs` の分子は **`Found/` にある条なしの `.src`** である
(`item := "Kind N.M"` **完全一致**)。★★**`"Definition 1.5, (i)"` のような
条番号つきは分子に入らない**——サブ項目もまた「条つき」である。

### ★★§1 の残作業(条ごと)

| 命題 | `Found/` に条なしで在る条 | 残る条 | 残りを埋めるのに要るもの |
|---|---|---|---|
| Definition 1.1 | — | 全体 | 射影埋め込み(`ArcModel` の構成) |
| Definition 1.2 | (ii) | (i) | 可逆層のテンソル積 or 移動補題 |
| Example 1.3 | **(i)** ←本日 | (ii) | コンパクト領域(非アルキメデス側) |
| Proposition 1.4 | (iv) | (i)(ii)(iii) | (i)(iii) 可逆層、(ii) 射影埋め込み |
| Definition 1.5 | (i)(iii)(iv) | **(ii)** | **Auslander–Buchsbaum** |
| Proposition 1.6 | — | 全体 | 射影埋め込み(`ArcModel` の構成) |

### ★★★誤り 1: 「`ℙⁿ` の点の関手で 5/24 → 7/24」は**過大**である

§9-18 と `memory/genell-track-b.md` にそう書いた。★**誤りである。**

`ℙⁿ` の点の関手が解くのは `Proposition 1.4, (ii)` と `Proposition 1.6` だが、
★★`Proposition 1.4` は **(i)(ii)(iii)(iv) が全部揃わないと数に入らない**。
(i)(iii) は可逆層の側が要るので、`ℙⁿ` だけでは埋まらない。

★★★**正しくは 5/24 → 6/24**(`Proposition 1.6` の 1 件のみ)。

### ★★★誤り 2: 「Auslander–Buchsbaum 1 本」は**規模を過小評価**していた

原文 `Definition 1.5, (ii)` の逐語:
> (ii) Let E ⊆ Z be an effective Cartier divisor contained in the regular locus of
> a normal noetherian scheme Z. Then observe that the closed subscheme Ered ⊆ Z
> is also an effective Cartier divisor.

★機構は「正則局所環 `𝒪_{Z,x}` は UFD だから `√(f)` が単項」である。
★★mathlib の在庫は **`RingTheory/RegularLocalRing/Defs.lean` 1 ファイルだけ**で、
Auslander–Buchsbaum は入っていない。

★★★**「1 本」だが、その 1 本が Serre の特徴づけ(正則 ⟺ 有限大域次元)
＋ Nagata の補題 ＋ 次元の帰納**である。**mathlib 級の企画**であって、
セッション単位の作業ではない。

★ただし**位置は良い**——`Definition 1.5` は (i)(iii)(iv) が既に揃っているので、
**この 1 本だけで 1 件動く**。§1 で単一の障害しか残っていないのはこれだけである。

### ★§3・§4 の残りも測った(初めて)

| 節 | 済 | 残り | 何が要るか |
|---|---|---|---|
| §3 | `Lemma 3.1` / `Lemma 3.6` | 7 件 | Faltings 高さ・半安定還元・`M_ell` の理論 |
| §4 | `Lemma 4.1` / `Remark 4.1.1` / `Lemma 4.2` | `Cor 4.3` / `Cor 4.4` | **論文の最終結果**(IUT 依存) |

★★§4 が 3/5 と見えるのは残り 2 件が**最終結果**だからで、近いのではない。
★`Lemma 3.5` / `Lemma 3.7` を「補題だから自己完結だろう」と当たったが、
**`htFalt` / `deg∞` / `M_ell` が要る**——mathlib に無い。

### ★★★結論: 24 件のうち安い項目は残っていない

**19 件すべてが、mathlib 級の基盤(Auslander–Buchsbaum / Faltings 高さ /
楕円曲線のモデュライ / `ℙⁿ` の点の関手 / 可逆層のテンソル積)か、
IUT そのものを要求する。**

★これは「進まない」の言い換えではなく、**どれを取るかの判断材料**である。
★★★単一障害で 1 件動くのは **`Definition 1.5`(Auslander–Buchsbaum)**と
**`Proposition 1.6`(射影埋め込み)**の 2 つだけである。

### ★本日の実装: `Example 1.3, (i)`(`Found/GenEll/DegSubset.lean`)

原文の (i) は**定義のみ**である。4 つの記号を全部作った:

| 記号 | Lean |
|---|---|
| `X(ℚ̄)^{≤d}` (`d ∈ ℕ ∪ {∞}`) | `leDeg X D (d : ℕ∞)` |
| `X(ℚ̄)^{=d}` | `eqDeg X D d` |
| `E^{≤d}` / `E^{=d}` | `subLeDeg E d` / `subEqDeg E d` |
| Galois-finite | `GaloisFinite E` |

★原文が付けている但し書き 2 つも定理にした:
- `X(ℚ̄)^{≤∞} = X(ℚ̄)` → `leDeg_top`
- 「最小定義体を見よ」→ `mem_eqDeg_iff`
  (`x ∈ X^{=d}` ⟺ **定義体の次数の最小値が `d`**)

★★★**商 `UPoint` の上で「次数」は代表元ごとに違う**——だから `≤ d` は
**「そういう代表元が存在する」**と読むのが正しい。存在量化なので
`Quot` の上で well-defined であることは自動である。

★これは `Lemma 3.7`(Galois-finite な例外集合 `Exc`)が直接要求するものである。

---

## §9-20 ★★★Arakelov 理論・Galois 表現のスケルトン——**個数を数えられるようにした**

**2026-08-17。**ユーザーの問い:
> S1からS4で現在、0/3，0/9，0/9，0/3という進捗ですか?
> これらを前にArakelov, Galois表現が存在しているのですか?
> Arakelov, Galois表現での形式化すべき個数は把握(スケルトン作成)できていますか?

★★★**3 つ目への答えは「できていなかった」である。**
`Interface/GenEll/` の `waiting` は**大きな文字列 1 本**で待っており、
「何本作れば埋まるか」が数えられなかった。★本節はそれを割った記録である。

### ★★★Arakelov 理論 —— **9 件**

`Interface/Arakelov/{LineBundle,ArcSpace,APic}.lean`

| # | obligation | 我々の到達点 | mathlib(2026-08-17 実測) |
|---|---|---|---|
| B1 | `Pic(X)`(可逆層の群) | — | ★`IsLocallyFree` は**ある**。**階数**と**層のテンソル積**が無い(前層版は有る) |
| B2 | `𝒪_X(D)`(Cartier → 可逆層) | ★`comap_mul` 証明済 | `IdealSheafData` / `comap` は有る |
| B3 | `Pic(Spec 𝓞_F) ≅ ClassGroup` | — | `ClassGroup` は有る。橋が無い |
| C1 | `X^arc` の位相・`ι_X`・コンパクト性 | ★★**構成済**(`ArcModel`) | ——(GAGA は**不要**と判明) |
| C2 | ℤ-固有 ⇒ 射影埋め込み | — | ★★★**層 C の律速**。`ℙⁿ` の点の関手が無い |
| C3 | 解析化 `L^arc` と hermitian 計量 | ★`IsConjInvariant` 定式化済 | (B1) に従属 |
| D1 | `APic(X)` | ★因子表示で対応物あり | (B1)+(C3) に従属 |
| D2 | `APic(Spec 𝓞_F)` と `deg_F` | ★★**`ADiv`/`deg_F`/`APrc` 実装済**、底変換不変性も証明済 | 橋が無い |
| D3 | `ht_L̄ : X(ℚ̄) → ℝ` と `Prop 1.4` | ★★★**`U_X(ℚ̄)` では構成済** | (B2) が入れば全域化される |

★★★**律速は B1(層のテンソル積)と C2(`ℙⁿ` の点の関手)の 2 本である。**

### ★★★Galois 表現 —— **8 件**

`Interface/GaloisRep/{Torsion,Representation,Reduction}.lean`

| # | obligation | mathlib / FLT(2026-08-17 実測) |
|---|---|---|
| G1 | `E[n] ≅ (ℤ/n)²` | ★★★**両方に無い**(FLT は sorry)。**S3 の最初の壁** |
| G2 | Tate 加群 `T_l E ≅ ℤ_l²` | 0 件。(G1) に従属し、入れば機械的 |
| G3 | `ρ_{E,l} : Gal → GL₂(ℤ_l)` | 行き先(`GeneralLinearGroup`)と定義域(`AlgEquiv`)は**書ける**。★**Weil 対**が無い |
| G4 | `mod l` 表現 | `PadicInt.toZMod` は有る。★`Lemma 3.1` は**我々が実装済** |
| G5 | 像が `SL₂` を含む(`Theorem 3.8`) | S3 の最後の段 |
| G6 | Tate 曲線と局所高さ(`Definition 3.3`) | ★★両方に無い。典拠 [FC] III, Cor 7.3 |
| G7 | 半安定還元と `𝓞_L` 上のモデル | `Reduction.lean` は有るが Néron モデルは無い |
| G8 | Faltings 高さ `ht^Falt = deg(ω_E)` | ★★★**Arakelov 側との合流点**((D2) が要る) |

★★★**`E[n] ≅ (ℤ/n)²` が入口である**——これが無いと `GL₂` の **`2`** が書けない。

### ★★捩れ部分群は posit しなかった

`E[n]` そのものは mathlib の `W.toAffine.Point`(`AddCommGroup`、`[DecidableEq K]` が要る)から
**今すぐ書ける**ので、`torsionPoints` は **`def` として定義**した。
★posit したのは**構造定理の方だけ**である。★★**posit を最小にする**という規律の適用である。

### ★★★数え方

`node tools/check.mjs` の「Interface 実装待ち」の行が**そのまま件数**になる。
2026-08-17 時点で **26 件**(Arakelov 9 + Galois 表現 8 + 既存 9)。

★★**増えたように見えるが、これは後退ではない**——
これまで 1 本の文字列に畳まれていたものが、**個別に埋められる形で数えられるようになった**。
★1 件埋めるごとにキューが減る。

### ★器具の穴を 3 つ見つけた

| # | 症状 | 原因 |
|---|---|---|
| 1 | `lake env lean` は通るのに `lake build` が落ちる | ★**auto-bound universe の設定が違う**。`universe u v` を明示して解決 |
| 2 | G2 が `PulledBackClassData..nonvacuous` を探す | ★宣言名の正規表現が `.{u, v}` を名前の一部として拾う |
| 3 | 逐語照合が 67/69 文字で止まる | ★引用を**字下げして**フィールドの docstring に置いた |

★★★**1 は規律への影響が大きい**——「`lake env lean` を通ったら即コミット」だけでは
`lake build` の失敗を取り逃す。★**コミット後に `lake build` を回し、落ちたら次のコミットで直す。**

---

## §9-21 ★★★原典を集めた——`Definition 1.5, (ii)` の見積りが変わった

**2026-08-17。**ユーザーの問い「必要な書籍はすべて有りますか?なければ公開記事や
解説サイトなどを調査して必要な情報を取得できますか?」に答えて集めた記録である。

### ★★★測定の更新: Auslander–Buchsbaum は「読む物が無い」段ではなかった

§9-20 で「Auslander–Buchsbaum は Serre の特徴づけ + Nagata の補題 + 次元の帰納で
**mathlib 級の企画**」と書いた。規模の判断は変えないが、★**前提が 2 つ変わった**:

1. ★★★**証明が手元にある**——`Stacks Project` tag **`0AG0`**(Lemma 15.123.2)
   「A regular local ring is a UFD」の**完全な証明**。次元の帰納で、依存も全部辿れる:

   | 依存 | 内容 |
   |---|---|
   | Algebra 10.106.2 | 正則局所環は整域 |
   | Algebra 10.106.3 | `R/(x)` が正則 |
   | Algebra 10.120.6 | Noether 整域が UFD ⟺ 高さ 1 の素イデアルが全部単項 |
   | Algebra 10.110.6 | ★**局所化が正則**(Serre の判定法。ここが最難) |
   | Algebra 10.78.2 | 有限表示 + 各点で階数 1 自由 ⇒ 可逆 |
   | More on Algebra 15.123.1 | |

2. ★★**最難の入力が既に Lean4 で形式化されている**——
   Guan–Hu, *Formalization of Auslander–Buchsbaum–Serre Criterion in Lean4*
   (arXiv 2510.24818、2025-10 投稿 / 2025-12 改訂)。
   Rees の定理・Auslander–Buchsbaum 公式・Ischebeck の定理・Cohen–Macaulay 加群・
   Hilbert の syzygy 定理を含む。
   ★**ただし「正則局所 ⇒ UFD」そのものは含まないと本文が述べている。**
   ★公開先(GitHub)は未確認であり、mathlib に入っているかも未確認。

★★★**したがって `Definition 1.5, (ii)` は「誰も書いていない」ではなく
「他人が半分やっており、残り半分は証明が手元にある」段である。**

### ★★取得した原典(15 本。すべて公開分)

| tag | 出典 | 効く obligation |
|---|---|---|
| `Szp` | Szpiro, Degrés, intersections, hauteurs (Ast. 127) | B1・C3・D2。★`[GenEll]` 自身の参考文献 |
| `MBMetr` | Moret-Bailly, Métriques permises (Ast. 127) | ★★C3・D1(エルミート計量) |
| `Elkik` | Elkik, Fonctions de Green, volumes de Faltings (Ast. 127) | ★★C3・G8(Green 関数) |
| `MBComp` | Moret-Bailly, Compactifications, hauteurs et finitude (Ast. 127) | D2・D3・`Prop 1.4 (iv)` |
| `RayHt` | Raynaud, Hauteurs et isogénies (Ast. 127) | ★★★G8。`Lemma 3.5` の機構そのもの |
| `DelRep` | Deligne, Représentations ℓ-adiques (Ast. 127) | ★★G3・G4 |
| `DelSB616` | Deligne, Bourbaki 616(Faltings の解説) | G7・G8。`[Falt]`(有料)の代替 |
| `RaySB427` | Raynaud, Bourbaki 427(Mumford 構成) | ★★★G6。`[FC]`(書籍)の代替 |
| `EGA1` | Grothendieck, EGA I (Publ IHES 4) | B1 |
| `EGA2` | Grothendieck, EGA II (Publ IHES 8) | ★★★C2。`Proj` と射影射 |
| `GS` | Gillet–Soulé, Arithmetic intersection theory (Publ IHES 72) | ★★D1・D2・D3。算術版の移動補題 |
| `SouleAI` | Soulé, Arithmetic Intersection(講義録) | `GS` の入口 |
| `MilneAV` | Milne, Abelian Varieties(講義録) | ★★G1・G2・G3。`[Serre]` の代替になる範囲 |
| `Stacks` | The Stacks Project(全 7654 頁) | ★★★`Def 1.5 (ii)`・B1・B2・C2 |
| `ABSLean` | Guan–Hu, ABS 判定法の Lean4 形式化 | ★★★`Def 1.5 (ii)` の最難の入力 |

★★取得元は **Numdam / IHES / arXiv / jmilne.org / Stacks** —— すべて公開である。

### ★★★取得できないもの(著作権)

| 参考文献 | 種別 | 代替 |
|---|---|---|
| `[FC]` Faltings–Chai | Springer 書籍 | ★`RaySB427`(Mumford 構成)——**同値ではない** |
| `[Serre]` Abelian ℓ-adic Representations | Benjamin 書籍 | ★`MilneAV` + `DelRep`——**範囲が狭い** |
| `[Silv1]` `[Silv2]` | Springer 書籍の章 | `Szp` / `MBComp` / `GS` |
| `[Edw]` Riemann's Zeta Function | Academic Press 書籍 | ★不要(`PrimeNumberTheoremAnd` で代替済) |
| `[Falt]` / `[Merel]` / `[Elkies]` / `[vF]` | 有料誌 | `DelSB616` が `[Falt]` を解説 |

★★**購入または機関アクセスが要る。**海賊版は取らない。

### ★★★登記は「未目視」と明記した

全 15 本とも `notationRisk := "unmeasured"` / `verifiedPages := []`。
テキスト層は取れている(Stacks は 1983 万字)が、
★**どの記号が壊れるかは未測定**である。逐語を `.src` に使う前に目視すること。

### ★`1_Structured` は作っていない

**理由**: 構造化は逐語の正確さを前提にするが、**まだ 1 ページも目視していない**。
★`[GenEll]` では pdftotext が**行列の出現順を入れ替える**という実害が出ており、
1985 年のフランス語スキャンや 7654 頁の Stacks を無検証で構造化すると、
**壊れた逐語が `.src` に固定される**。
★★**実際に引く箇所だけを目視 → その節だけ構造化**、という順が正しい。

---

## §9-22 ★★★`mathlib4_fork` を実測した——依存には加えられないが、参照先になる

**2026-08-17。**arXiv 2510.24818 の PDF から GitHub の所在を抽出し
(本文はハイパーリンクなので PDF の圧縮ストリームから URL を復元した)、
**clone せず raw ファイルを取って計数**した。

**`https://github.com/Thmoas-Guan/mathlib4_fork`**(公開・Apache-2.0・2026-08-03 更新)

### ★★★得られたもの: `isRegularLocalRing_localization` が sorry 無しで在る

ブランチ `ABS-Criterion-Project-new` の
`Mathlib/RingTheory/RegularLocalRing/` に 5 ファイル・**1313 行・sorry 0**:

| ファイル | 行 | 中身 |
|---|---|---|
| `Defs.lean` | 121 | 手元 mathlib にもあるが拡張されている |
| `Basic.lean` | 329 | ★手元 mathlib に**無い** |
| `GlobalDimension.lean` | 122 | ★手元 mathlib に**無い** |
| `AuslanderBuchsbaumSerre.lean` | 694 | `IsRegularLocalRing.of_globalDimension_lt_top` / `generate_by_regular` |
| `Localization.lean` | 47 | ★★★`Auslander_Buchsbaum_Serre` / **`isRegularLocalRing_localization`** |

★★★**`isRegularLocalRing_localization` は Stacks 10.110.6 そのもの**であり、
§9-21 で「`Definition 1.5, (ii)` の最難の依存」と測ったものである。**それが sorry 無しで在る。**

### ★★★依存には加えられない(技術的に不可能)

**これは mathlib の*フォーク*である。**mathlib と同時に `require` すると
`Mathlib.*` のモジュール名が衝突する。★`PrimeNumberTheoremAnd` のように
「別ライブラリとして足す」ことはできない。

★★したがって道は **(a) 参照して自分で書く** か **(b) 該当ファイルを移植する** の 2 つ。

### ★移植する場合の実測

- toolchain: fork は **`v4.26.0-rc2`**、我々は **`v4.31.0`** —— ★**5 版の drift** を吸収する必要がある。
- ライセンス: Apache-2.0。★移植するなら `Copyright (c) 2025 Nailin Guan` の表示を残すこと。
- 最短は `Localization.lean`(47 行)だが、`AuslanderBuchsbaumSerre.lean`(694 行)に依存する。

### ★★★それでも `Definition 1.5, (ii)` は埋まらない

186 本のブランチ名を全部見たが、**UFD / factorization 系は 0 件**。
論文本文も「正則局所 ⇒ UFD は含まない」と述べている。★残る 3 本は依然として無い:

| Stacks | 内容 | mathlib | fork |
|---|---|---|---|
| 10.110.6 | 局所化が正則 | ✗ | ★★★**在る(sorry 0)** |
| 10.106.2 | 正則局所環は整域 | ✗ | ✗ |
| 10.120.6 | Noether 整域: UFD ⟺ 高さ1素が単項 | ✗ | ✗ |
| 15.123.2 | **正則局所環は UFD** | ✗ | ✗ |

★★★**「最難の 1 本が他人によって済んでいる」段であって、「終わっている」段ではない。**

### ★保存した場所

`ResearchPaper/0_Source/mathlib4_fork-ABS/`(5 ファイル + `LICENSE` + `lean-toolchain`)。
実測は `ResearchPaper/lean-ecosystem.json` の
`mathlib4_fork (ABS-Criterion-Project)` に登録した。

---

## §9-23 ★★★C1 の実装 —— 穴は 1 点まで絞れた(次はここから)

**2026-08-17。**退化封じ(§9-22 の後)に続き、**実際に埋める作業**として C1
(`Interface/Arakelov/ArcSpace.lean` の `ArcSpaceData`)を実装した。
`Found/Arakelov/` に **9 ファイル・すべて sorry 0**。

### ★★C1 の 7 要求の状態

| 場 | 実装 | 場所 |
|---|---|---|
| `Arc X` | ✅ | `complexPoints`(既存) |
| `equivComplexPoints` | ✅ | `Equiv.refl` |
| `evalAffine` | ✅ | `ArcEval.lean`(`Spec.preimage`) |
| `conj` / `conj_involutive` | ✅ | `ArcConjInvol.lean` |
| `conj_continuous` | ✅ | `ArcTopology.lean` |
| `topology` | ✅ | `ArcTopology.lean`(`⨆` over affine opens) |
| `topology_affine` | ✅ | `ArcTopologyAffineEq.lean` |
| `topology_openImmersion` | ★**前半のみ** | `ArcOpenImmersion.lean` |

### ★★★残る穴は 1 点である

    induced (· ≫ f) (arcTopology Y) ≤ arcTopology X

すなわち「**`(· ≫ f)` が開写像**」。★手元にある部品:

| 部品 | 場所 |
|---|---|
| `continuous_comp_openImmersion`(前半) | `ArcOpenImmersion.lean` |
| `comp_openImmersion_injective` | 同 |
| `preimage_image_comp_openImmersion` | 同 |
| `isOpen_arcBasicOpen`(`{p ǀ g(p) ≠ 0}` が開) | `ArcBasicOpen.lean` |
| `imageAffineOpen` と mathlib の `affineOpensEquiv` の一致 | `ArcOpenImmersion.lean` |

### ★★★★次のセッションの 3 段(ここから始める)

1. **`W := (· ≫ f) '' V` が `arcTopology Y` で開**であることを示す。
   `arcTopology Y = ⨆` なので、各アフィン開 `V' ⊆ Y` について
   `chart_{V'}⁻¹ W` が `arcTopologyOpen V'` で開であることに帰着する。
2. `chart_{V'}⁻¹ W = {r : Arc V' ǀ r ≫ V'.ι が f を経由し、その因子が V に入る}`。
   ★★**「`f` を経由する」は「像の点が `f.opensRange` に入る」**であり、
   `V' ∩ f.opensRange` は `V'` の開集合。
   ★アフィン `V'` では基本開集合の合併なので **`isOpen_arcBasicOpen` で開**。
3. その開集合の上で**逆写像 `ψ : {r ǀ …} → Arc X`** を作り
   (`r ≫ V'.ι = ψ(r) ≫ f`、一意性は `comp_openImmersion_injective`)、
   **`ψ` の連続性**を示す。★★`chart_{V'}⁻¹ W = ψ⁻¹ V` となって終わる。

★★**鍵は `IsOpenImmersion.affineOpensEquiv`**(mathlib、順序同型)——
`X` のアフィン開被覆と `Y` の被覆のうち像に入るものが 1 対 1 に対応する。

### ★実装で 2 度かかった罠(記録)

★★★**`simp only [Function.comp_def]` は `Scheme.Hom.mk` まで展開して壊す。**
`rw [h]` を使い、**`h` を「合成の形」で述べる**(ゴールが既に合成形なので)。

★位相のインスタンスは `letI` で**明示的に**入れる(`⨆` の各成分は自動で決まらない)。

★★★**`TopologicalSpace` の `≤` は「細かい」**。`le_def` の表示は対称で読めないので、
**`⊥` が離散かを試して**確定させること(2026-08-17 に一度誤った)。

---

## §9-24 ★★★C1 の到達点(2026-08-17 終)—— 残るのは 1 点、その論法まで特定した

`Found/Arakelov/` は **11 ファイル・すべて sorry 0**。

### ★★C1(`ArcSpaceData`)の 7 要求

| 場 | 状態 | 場所 |
|---|---|---|
| `Arc X` / `equivComplexPoints` | ✅ | `complexPoints`(既存)/ `Equiv.refl` |
| `evalAffine` | ✅ | `ArcEval.lean` |
| `conj` / `conj_involutive` | ✅ | `ArcConjInvol.lean` |
| `conj_continuous` | ✅ | `ArcTopology.lean` |
| `topology` | ✅ | `ArcTopology.lean`(`⨆` over affine opens) |
| `topology_affine` | ✅ | `ArcTopologyAffineEq.lean` |
| `topology_openImmersion` | ★**1 点残** | `ArcOpenImmersion.lean` / `ArcLift.lean` |

### ★★★`topology_openImmersion` の内訳

    topology X = induced (· ≫ f) (topology Y)      (f : X ⟶ Y が開埋め込み)

| 向き | 状態 | 定理 |
|---|---|---|
| `topology X ≤ induced` | ✅ | `continuous_comp_openImmersion` |
| `induced ≤ topology X` | ★**1 点残** | 下記 |

残り 1 向きは「**`(· ≫ f)` が開写像**」。★手元の部品:

| 部品 | 定理 |
|---|---|
| chart が chart に収まる | `imageAffineOpenIso_fac` |
| アフィン開の 1 対 1 対応(mathlib の順序同型) | `imageAffineOpen_eq_affineOpensEquiv` / `exists_imageAffineOpen` |
| 単射性 | `comp_openImmersion_injective` |
| 基本開集合は開 | `isOpen_arcBasicOpen` |
| **点が `D(g)` に落ちる ⟺ `g(r) ≠ 0`** | `mem_basicOpen_base_iff` |
| **像の点を取る写像は連続** | `continuous_base_default` |
| **点が開集合に落ちる条件は開** | `isOpen_landsIn` |
| **像への所属 ⟺ 点が像に落ちる** | `mem_range_comp_iff_base_default` |
| 持ち上げの一意性 | `lift_eq_of_comp_eq` |

### ★★★★残る 1 点と、その論法(次のセッションはここから)

示すべきは、`V` が `arcTopology X` で開のとき

    ∀ V' ∈ Y.affineOpens, chart^Y_{V'} ⁻¹' ((· ≫ f) '' V)  が arcTopologyOpen V' で開

★★**難所は「`f⁻¹ᵁ V'` がアフィンでない」ことである。**
`arcTopology` はアフィン chart で定義したので、
`Arc (f⁻¹ᵁ V') → Arc X` の連続性が直接には使えない。

★★★**論法(アフィン細分)**:

1. `W' := V' ⊓ f.opensRange`(`V'` の開)。像に入る `r` は `W'` に落ちる
   ——`mem_range_comp_iff_base_default` + `isOpen_landsIn` で**開集合**である。
2. `W'` を **`Y` のアフィン開 `V''` で被覆**する(`V'' ≤ W'`)。
   ★`exists_imageAffineOpen` により各 `V''` は `X` のアフィン開 `U''` の像である。
3. `{r ∈ Arc V' ǀ r が V'' に落ちる}` の上では
   `chart^Y_{V'}⁻¹((· ≫ f)''V) = (chart^X_{U''} を経た像)` になり、
   `chart^X_{U''}⁻¹ V` が開(`V` が開だから)であることに帰着する。
4. 開集合の**合併**なので全体も開。

★★見積り: **100〜200 行**。★部品は 9 つとも揃っているので、
新たな mathlib 探索は要らない見込みである。

### ★★実装で 3 度かかった罠(必ず読むこと)

1. ★★★**`simp only [Function.comp_def]` は `Scheme.Hom.mk` まで展開して壊す。**
   `rw [h]` を使い、**`h` を「合成の形」で述べる**(ゴールが既に合成形)。
2. ★★位相のインスタンスは `letI` で**明示的に**入れる(`⨆` の成分は自動で決まらない)。
3. ★★★**`TopologicalSpace` の `≤` は「細かい」**。`le_def` の表示は対称なので、
   **`⊥` が離散かを試して**確定させる(一度誤って `⨅`/`⨆` を逆にした)。
4. ★`Spec ℂ` の点の型は **`↥(Spec (CommRingCat.of ℂ))`** で
   `PrimeSpectrum (…)` と**構文的に別**(defeq だが `rw` は噛まない)。

### ★★★★論法の底まで詰めた(2026-08-17 追記)—— 基底事例は「局所化の延長の連続性」

§9-24 の 4 段を細かく追うと、**持ち上げの連続性は再帰的に現れる**が、
★★アフィン細分をすると**アフィン間の持ち上げ**に落ちる。そこで底に着く。

★具体的に: `V' = Spec B` がアフィン、`V'' = D(g) ⊆ V'` が基本開集合のとき

    Arc (Spec B) = Hom_Ring(B, ℂ)      (`ArcEval.lean` の `evalHom_injective` / `evalHom_Spec_map`)
    Arc (D(g))   = Hom_Ring(B_g, ℂ)

であり、持ち上げは **`φ : B → ℂ`(`φ(g) ≠ 0`)を `B_g → ℂ` へ延長する操作**である。

★★★**したがって残る 1 点の底は次の 1 文である**:

> `{φ ∣ φ(g) ≠ 0}` の上で、局所化の延長 `B_g → ℂ` は `φ` について連続である。

★機構は `b/gⁿ ↦ φ(b)/φ(g)ⁿ` であり、**分母が消えない所での除算の連続性**
(`ContinuousAt.div`)に帰着する。★mathlib の `IsLocalization.lift` が延長を与える。

★★これで**新たな数学は要らない**ことが確定した——残るのは
(a) 局所化の延長の連続性(上記 1 文)、(b) それを開被覆で貼る段、の 2 つである。

★★★見積りの更新: **100〜200 行**のうち、(a) が 40〜60 行、(b) が 60〜140 行。

---

## §9-25 ★★★★C1 の基底事例が閉じた —— 一般の場合の手順を精密化した

`Found/Arakelov/` は **12 ファイル・全ファイル sorry 0**。

### ★★★★取れた基底事例

    arcTopologyAffine (Localization.Away s)
      = induced (· ≫ awayι A s) (arcTopologyAffine A)

★★★これは `topology_openImmersion` の**アフィン(基本開集合)の場合そのもの**である。
★4 条件が揃った帰結:

| 条件 | 定理 |
|---|---|
| `awayLift → 合成` でもとに戻る | `awayLift_comp_awayι` |
| `合成 → awayLift` でもとに戻る | `awayLift_comp_awayι_self` |
| `awayLift` は連続 | `continuous_awayLift` |
| 合成は連続 | `continuous_comp_affine` |

★★**両向き連続な全単射なので部分空間位相になる。**

### ★★★★一般の場合の手順(次のセッションはここから、機械的に実行できる)

示すべきは、`V` が `arcTopology X` で開のとき `W := (· ≫ f) '' V` が
`arcTopology Y` で開であること。`isOpen_iSup_iff` で `V' ∈ Y.affineOpens` ごとに落とす。

1. **`V' ⊓ f.opensRange` を基本開集合に分解する。**
   `V'` はアフィンなので `isoSpec` で `Spec Γ(Y,V')` に移し、
   `PrimeSpectrum.isTopologicalBasis_basic_opens` で `⋃ D(gᵢ)` と書く。
2. **各 `D(gᵢ)` は `X` のアフィン開の像である。**
   `D(gᵢ) ≤ f.opensRange` なので `exists_imageAffineOpen` により
   `∃ Uᵢ ∈ X.affineOpens, f ''ᵁ Uᵢ = D(gᵢ)`。
3. **`{r ∈ Arc V' ǀ r が D(gᵢ) に落ちる}` の上で基底事例を当てる。**
   ★この集合は `isOpen_landsIn`(アフィン、`isoSpec` で輸送)により**開**。
   ★★その上では `arcTopologyAffine (Away gᵢ) = induced …`(基底事例)により
   位相が一致するので、`chart_{V'}⁻¹ W` は
   `chart^X_{Uᵢ}⁻¹ V`(`V` が開だから開)の像に一致する。
4. **合併を取る。** `chart_{V'}⁻¹ W = ⋃ᵢ (上の開集合)` なので開。

★★★**新たな数学は無く、既存 14 定理の組み合わせである。**
見積り **150〜250 行**、うち手順 1 の分解と手順 3 の輸送が主な作業。

### ★★使える定理の一覧(すべて `Found/Arakelov/`、sorry 0)

`evalHom_injective` / `evalAffine_comp` / `continuous_comp_affine` /
`continuous_evalAffine` / `t2Space_arcAffine` / `conjPoint_conjPoint` /
`evalAffine_conjPoint` / `continuous_conjPoint` / `arcTopology_spec` /
`continuous_chart_affine` / `imageAffineOpenIso_fac` / `exists_imageAffineOpen` /
`comp_openImmersion_injective` / `continuous_comp_openImmersion` /
`isOpen_arcBasicOpen` / `continuousOn_div_pow_evalAffine` /
`mem_basicOpen_base_iff` / `continuous_base_default` / `isOpen_landsIn` /
`mem_range_comp_iff_base_default` / `continuous_base_default_scheme` /
`isOpen_landsIn_scheme` / `awayLiftHom_mk` / `continuous_awayLift` /
`awayLift_comp_awayι` / `awayLift_comp_awayι_self` / `arcTopologyAffine_away`


---

## §9-26 (2026-08-17) ★★★★★C1 達成と、残り 16 件の**到達可能性の実測**

### ★★★★★C1 は落ちた —— Arakelov **1/9**

`ArcSpaceData.nonvacuous` が取れた(`Found/Arakelov/` 14 ファイル、sorry 0、標準 3 公理)。
`check.mjs` の「Interface 実装待ち」は **26 → 25 件**。

3 段の積み上げ:

| 段 | 定理 | 内容 |
|---|---|---|
| A | `arcTopology_opens_of_affine` | `Spec A` の**任意の**開部分スキーム |
| B 前半 | `isOpenMap_comp_of_isAffine` | 任意のアフィン標的 |
| B 後半 | `arcTopology_openImmersion` | 一般の開埋め込み |

★★★段 B の核心は「**`U ⊓ O` を経由する**」ことだった:

    (· ≫ U.ι) ⁻¹' ((· ≫ O.ι) '' V) = (· ≫ homOfLE) '' ((· ≫ homOfLE) ⁻¹' V)

`arcTopology = ⨆`(アフィン chart)なので開性は chart ごとに降ろせ、
各 chart は**アフィン**だから段 A が効く。

★★★**GAGA は要らなかった**——当初の見積りは誤りで、実際に要ったのは
**商位相と多項式の連続性**だけだった。

### ★★★★退化封じの見落とし(2 度目の教訓)

C1 は当初、次を**通していた**:

    evalAffine := fun _ _ _ => 0 ;  topology := fun _ => ⊤(密着)

`induced (fun _ _ => 0) Pi = ⊤` なので `topology_affine` は「`topology = ⊤`」を
要求するだけになり、`topology_openImmersion` も `induced g ⊤ = ⊤` で自明に成り立つ。

★`equivComplexPoints` は**台の型**を固定するが、**`evalAffine` の値**は固定しない。
★`evalAffine_spec`(評価は `Spec.preimage` が与える環準同型)を足して塞いだ。
★負の対照は `Check/Arakelov/ArcSpaceNondegenerate.lean`(離散・密着の両方が落ちる)。

★★★**教訓: 退化封じでは「値が自由に選べるフィールド」を列挙する。型の固定だけでは足りない。**

### ★★★★★残り 16 件は**すべて mathlib 規模の基礎**に閉ざされている(実測)

2026-08-17 に mathlib を全走査して測った結果:

| 障壁 | 実測 | 閉ざす obligation |
|---|---|---|
| **層の圏のモノイダル構造** | `SheafOfModules` に `MonoidalCategory` インスタンス **0 件** | B1 → B2 B3 C3 D1 D3(**6 件**) |
| **ℙⁿ の点の関手** | `ProjectiveSpectrum` に `Hom(Spec ℂ, Proj) ≅ ℙ(ℂ^{n+1})` **無し** | C2 |
| **E[n] ≅ (ℤ/n)²** | mathlib 無し、**FLT も sorry** | G1 → G2 G3 G4 G5(**5 件**) |
| **Tate 曲線** | mathlib 無し、FLT は 20 行の入口だけ | G6 → G8 |
| **Néron モデル** | mathlib 無し | G7 |

★★★**したがって「9/9・8/8」は、現在の mathlib からは 1 セッションでは到達しない。**
16 件はいずれも「我々が書き落とせば済む定理」ではなく、**基礎理論の建設**である。

### ★★★★B1(Pic)への道筋 —— **部品は在る**

B1 は 6 件を解くので投資先として最も効く。★2026-08-17 の実測で、
**必要な機械は mathlib に揃っていた**(組み立てられていないだけ):

| 部品 | 場所 | 状態 |
|---|---|---|
| 前層の対称モノイダル構造 | `ModuleCat/Presheaf/Monoidal.lean` | ★`monoidalCategory` / `symmetricCategory` **在り** |
| 層化(前層 → 層) | `ModuleCat/Presheaf/Sheafification.lean` | ★`sheafificationAdjunction`、**counit は iso**(反射的) |
| 局所化のモノイダル移送 | `CategoryTheory/Localization/Monoidal/` | ★`LocalizedMonoidal` **在り**(`W.IsMonoidal` が要る) |
| モノイダル圏の Picard 群 | `CategoryTheory/Monoidal/Skeleton.lean` | ★`Skeleton` の `CommMonoid`、その `ˣ` が群 |
| 局所自由性 | `ModuleCat/Sheaf/LocallyFree.lean` | ★`IsLocallyFree` **在り** |
| アフィンでの対応 `M^~` | `AlgebraicGeometry/Modules/Tilde.lean` | ★★`tilde` / `SpecModulesToSheafFullyFaithful` / `fromTildeΓ` **在り** |
| 環の Picard 群 | `RingTheory/PicardGroup.lean` | ★`CommRing.Pic` / `ClassGroup.equivPic` **在り** |

#### ★2 つの道

**道 1(一般):** 局所同型のなす `W` が `IsMonoidal` であることを示し、
`LocalizedMonoidal` で層の圏にモノイダル構造を移す。★正攻法だが重い。

**道 2(近道):** ★★★**可逆層に限れば層化は要らない**——
局所自由有限階数の層どうしの**前層テンソル積は既に層である**
(局所的に `𝒪^n` の有限直和で、層であることは局所的性質だから)。
★これなら `Localization/Monoidal` を通さずに済む。

#### ★残る難所

いずれの道でも `equivPicRing : Pic (Spec R) ≃* CommRing.Pic R` が要る。
★`Tilde.lean` の `fromTildeΓ` が**準連接層で iso** であることを言えば出る
——mathlib には `Quasicoherent.lean` が在るが、**同値そのものは無い**。

### ★★次の 1 手

★★★**道 2 の第 1 ブロック**:「局所自由階数 1 の層 2 つの前層テンソル積は層である」。
★これが取れれば `Pic X` の群構造が書け、B1 の 6 件が動き出す。


---

## §9-28 (2026-08-17) ★★★★B1 の本丸は **内部 Hom** だと確定した

### ★取れたもの(累計、いずれも sorry 0)

| `CommGroup (Pic X)` の公理 | 状態 | 与えるもの |
|---|---|---|
| 乗法 `⊗` | ✅ | `tensorModules`(前層でテンソル → 層化) |
| 可換 | ✅ | `tensorModulesComm`(前層の braiding を層化で送る) |
| 単位元 | ✅ | `tensorUnitLeft` / `tensorUnitRight`(`λ_` / `ρ_` + counit iso) |
| **結合律** | ★★★**残** | 内部 Hom が要る |
| **逆元** | ★★★**残** | 内部 Hom が要る(双対 `Hom(L, 𝒪)`) |

### ★★★★なぜ内部 Hom なのか(mathlib の証明を読んで判明)

結合律は `sheafify (sheafify (A ⊗ B).val ⊗ C) ≅ sheafify ((A ⊗ B) ⊗ C)` に帰着し、
これは単位 `η : (A ⊗ B) ⟶ sheafify(A ⊗ B).val` について

    η ▷ C ∈ W        (W = 層化が反転させる射のクラス)

を言うことである。★これは `MorphismProperty.IsMonoidal` の whiskering 条件そのもの。

★★★**mathlib は同じ命題を `Sites/Monoidal.lean` で証明しているが、
その道具は「モノイダル閉性」である**(2026-08-17 に証明本文を読んで確認):

    W.whiskerLeft : Hom(F ⊗ G, H) ≅ Hom(G, [F, H]) と
                    「H が層なら [F, H] も層」(isSheaf_functorEnrichedHom)

★★しかしそれは **値が固定のモノイダル圏 `A` に値をとる前層** `Cᵒᵖ ⥤ A` の話で、
**`PresheafOfModules`(係数が環の層)には適用できない**
——テンソルが `ℤ` 上でなく `R` 上だから。

### ★★★実測: 内部 Hom は mathlib に無い

    Algebra/Category/ModuleCat/Presheaf/  に Closed.lean 無し
    internalHom / MonoidalClosed / ihom / homObj  いずれも 0 件(Presheaf/ と Sheaf/ 全走査)

★逆元(双対 `L^∨ = Hom(L, 𝒪)`)にも同じものが要る。
★★★**したがって「`Hom_R(F, G)` の前層版とその随伴」の構築が B1 の本丸である。**

### ★次にやること(3 段)

1. `homObj F G : PresheafOfModules R₀` —— 切断は `Hom(F|_U, G|_U)`
2. 随伴 `Hom(F ⊗ G, H) ≅ Hom(G, homObj F H)`
3. 「`H` が層なら `homObj F H` も層」

★これが取れれば `W.IsMonoidal` → `LocalizedMonoidal` → `MonoidalCategory X.Modules`
→ `Pic X` の `CommGroup` → **B1 とそれに従属する 6 件**(B2 B3 C3 D1 D3)。

★★規模の見積り: mathlib の `Sites/Monoidal.lean` + `Enriched/FunctorCategory` に
相当する量。**mathlib PR 級**であり、1 セッションでは終わらない。

### ★★★Galois 側(G1-G8)の状況は変わらず

G1(`E[n] ≅ (ℤ/n)²`)が 5 件を閉ざしている。mathlib にも FLT にも無い(FLT は sorry)。
★`TorsionStructureData` は 2 条(`structure_eq` と `torsion_finite`)からなり、
**witness には両方が要る**ので、片方だけでは件数は動かない。


---

## §9-29 (2026-08-17) ★★★★★B1 —— 5 ブロック積んだ。残りは局所全単射の 2 条

### ★取れたもの(累計、すべて sorry 0)

| ブロック | ファイル | 内容 |
|---|---|---|
| 1 | `PicPresheafTensor.lean` | `MonoidalCategory` / `SymmetricCategory` on `X.PresheafOfModules` |
| 2 | `PicSheafTensor.lean` | `tensorModules`、`tensorModulesComm`(**可換**) |
| 3 | 同上 | `tensorUnitLeft` / `tensorUnitRight`(**単位元**) |
| 4 | `PicRestrict.lean` | `restrictOpen`(**開集合への制限**)、`restrictOpenIso` |
| 5 | `PicRestrictTensor.lean` | ★**制限はテンソル積と両立**(`Functor.Monoidal.μIso`) |

### ★★★★見積りを 3 度外した(記録)

| 何 | 当初の見積り | 実際 |
|---|---|---|
| C1 の `X^arc` | 複素解析空間と GAGA が要る | ★商位相と多項式の連続性だけ |
| 開集合への制限 | 自作せねばならない | ★`SheafOfModules.over` が在った |
| 制限とテンソルの両立 | 自作せねばならない | ★`pushforward₀OfCommRingCat.Monoidal` が在った |

★★★**「無い」と決める前に測る。**3 度効いた。

### ★★★★残るのは 2 条だけ

結合律は `W (η ▷ C)`(`η` は層化の単位)に帰着し、
`WEqualsLocallyBijective` により **局所全射 + 局所単射**の 2 条になる。

#### ★(a) 局所全射 —— **局所自由性は要らない**、直接証明できる

`x ∈ (Q ⊗ C)(U)` は純テンソルの有限和 `Σ qᵢ ⊗ cᵢ`。
各 `qᵢ` について `imageSieve η qᵢ ∈ J U` なので、有限個の交わりも被覆
(`J.intersection_covering`)。その上で `x` は像に入る。

★★**証明の型**: 「`imageSieve (η ▷ C) x ∈ J U` を満たす `x` の集合」が
部分加群であって純テンソルを含むことを言い、`TensorProduct.induction_on` で閉じる。
- `0` … `imageSieve` は `⊤`
- 和 … 2 つの被覆篩の交わり
- 純テンソル … `imageSieve_mem`

#### ★(b) 局所単射 —— **ここで局所自由性が要る**

一般の `C` では偽(平坦性が要る)。★可逆層は**局所的に `𝒪`** なので平坦であり、
第 4・第 5 ブロック(制限と、制限がテンソルと両立すること)を使って
局所自明化の開集合上へ落とす。

### ★★次の 1 手

★★★**(a) 局所全射から書く。**局所自由性を使わないので独立に取れる。
その後 (b) に `IsLocallyFree` の局所自明化を当てる。


---

## §9-30 (2026-08-17) ★★★★★★B1 —— 14 ブロック積み、`CommGroup` の 5 公理が揃った

### ★取れたもの(すべて sorry 0、`Found/Arakelov/` 7 ファイル)

| # | ファイル | 内容 |
|---|---|---|
| 1 | `PicPresheafTensor` | 前層加群の対称モノイダル構造(橋渡し) |
| 2–3 | `PicSheafTensor` | `tensorModules` / 可換 / 単位律 |
| 4 | `PicRestrict` | 開集合への制限 |
| 5 | `PicRestrictTensor` | 制限とテンソルの両立(mathlib の `Monoidal` 在庫) |
| 6 | `PicLocalSurj` | 局所全射性はテンソルで保たれる(★局所自由性不要) |
| 7–8 | `PicLocalInj` / `PicLocality` | 局所単射性の基底・輸送・**局所性**(mathlib に無かった) |
| 9 | `PicLocalBasis` | 局所基底から局所単射性(★`Over` 不要) |
| 10 | `PicWhiskerW` | `W (f ▷ M)` / `W (M ◁ f)` |
| 11–13 | `PicAssoc` | `W` → `IsIso` → 可逆層とのテンソル → **結合律** |
| 14 | `PicType` | 可逆層の型 / `one` / **`symm`(逆元)** |

### ★★★`CommGroup (Pic X)` の 5 公理

| 公理 | 定理 |
|---|---|
| 乗法 | `tensorModules` |
| 可換 | `tensorModulesComm` |
| 単位元 | `tensorUnitLeft` / `tensorUnitRight` / `InvertibleSheaf.one` |
| 結合律 | `tensorModulesAssoc` |
| 逆元 | `InvertibleSheaf.symm` |

### ★★★★残り 2 点(いずれも「局所自明化を presheaf の同型に上げる」が要る)

#### (a) テンソル積が可逆層で閉じること

`IsLocallyRankOne (tensorModules L M).val` が要る。★問題は
**層化された**テンソル積の切断が直接計算できないことである。

★★道筋: `L.val|_V ≅ 𝒪|_V` と `M.val|_V ≅ 𝒪|_V` を
**`Over V` 上の前層の同型**として持てば、第 5 ブロック(制限はモノイダル)で
`(L.val ⊗ M.val)|_V ≅ 𝒪|_V` となり、これは**層**なので層化が効かない。

★★★したがって `IsLocallyRankOne` を
「各 `V` で基底 1 本」から「`Over V` 上の前層同型」へ**強める**必要がある。
★第 9 ブロック(局所単射性)は弱い形で十分だったが、こちらは強い形が要る。

#### (b) 同型類への商と `CommGroup` インスタンス

`Pic X := Quotient (可逆層を同型で割ったもの)` とし、
乗法の well-defined 性(`tensorModulesIso`——第 2 ブロックに在る)と
5 公理を商へ降ろす。★★`Interface` は `Pic : Scheme → Type`(**小さい型**)を
要求するので `Shrink` が要る——本質的小ささの証明が別途要る。

### ★★★この turn で 5 度踏んだ罠(記録)

★★★**インスタンス束縛子は「型の書き方の違い」をまたげない。**

    X.PresheafOfModules   と   PresheafOfModules (X.presheaf ⋙ forget₂ _ _)

は `rfl` で等しいが、探索は片方でしか成功しない。
★**対処: 局所単射性・局所全射性は `[..]` ではなく `(h : ..)` で受ける。**
★もう 1 つ: **依存を書き換えたら `lake build` で先に olean を作り直す**
——さもないと名前付き引数が古い署名で解決される。


---

## §9-31 (2026-08-17) ★★★★★★★B1 —— 17 ブロック、`Pic X` が可換群になった

### ★取れたもの(`Found/Arakelov/` 10 ファイル、すべて sorry 0)

| # | ファイル | 内容 |
|---|---|---|
| 1 | `PicPresheafTensor` | 前層加群の対称モノイダル構造(橋渡し) |
| 2–3 | `PicSheafTensor` | 層加群のテンソル積・可換・単位律 |
| 4–5 | `PicRestrict` / `PicRestrictTensor` | 開集合への制限・制限とテンソルの両立 |
| 6 | `PicLocalSurj` | 局所全射性はテンソルで保たれる |
| 7–8 | `PicLocalInj` / `PicLocality` | 局所単射性の基底・輸送・局所性 |
| 9 | `PicLocalBasis` | 局所基底から局所単射性 |
| 10–13 | `PicWhiskerW` / `PicAssoc` | `W` → `IsIso` → **結合律** |
| 14 | `PicType` | 可逆層の型・逆元 |
| 15 | `PicLocalTrivial` | **局所自明性のテンソル閉性** |
| 16 | `PicGroup` | ★**可逆前層の群構造**(層化を通さない道) |
| 17 | `PicQuotient` | ★★**`PicPre X` と `CommGroup`** |

### ★★★★★層化を通さない道が決定打だった

第 13 ブロックまでは「前層でテンソル → 層化」で組んでいた。
★★★**層化を通すたびに「局所自明性が保たれるか」を示さねばならない**。

★★★★第 16 ブロックで**前層の段に移した**ところ:

| 何 | 前層の段 | 層化を通す段 |
|---|---|---|
| モノイダル構造 | ★mathlib が持っている | 無い(13 ブロックかけて作った) |
| 結合律・単位律・可換律 | ★★**無料**(`α_` / `λ_` / `ρ_` / `β_`) | 局所論法が要る |

★**局所自明な前層は自動的に層**(層条件は局所的)なので数学的に一致する。

### ★★★★B1 の残り 3 点

#### (1) 引き戻し `f^* : Pic Y → Pic X`

★mathlib は `Scheme.Modules.pullback f`(**層**の側)と
`pullbackId` / `pullbackComp` / `PullbackFree.lean` を持つ(2026-08-17 実測)。
★★★**しかしモノイダル性(`f^*(A ⊗ B) ≅ f^*A ⊗ f^*B`)は未登録**である。
★また我々の `Pic` は**前層**の側にあるので、まず層の側へ渡す必要がある。

#### (2) 局所自明な前層は層である

★`Interface` の `sheafOf` に繋ぐのに 1 回だけ要る。
「層条件は局所的」を site の言葉で言う必要がある。

#### (3) `equivPicRing` —— アフィンで `CommRing.Pic` と一致

★`AlgebraicGeometry/Modules/Tilde.lean` の `tilde` /
`SpecModulesToSheafFullyFaithful` / `fromTildeΓ` が材料。

### ★★★この turn の教訓(5 度踏んだ)

★★**インスタンス束縛子は「型の書き方の違い」をまたげない。**
`X.PresheafOfModules` と `PresheafOfModules (X.presheaf ⋙ forget₂ _ _)` は
`rfl` で等しいが探索は片方でしか成功しない。
★対処: 条件を `[..]` でなく `(h : ..)` で受ける。
★`Nonempty` から `obtain` して**データ**は返せない——結論も `Nonempty` にする。
★依存を書き換えたら `lake build` で先に olean を作り直す。


---

## §9-32 (2026-08-17) ★★★★★★B1 —— 18 ブロック、Interface と直結。残り 2 点の材料を実測した

### ★第 18 ブロック: Interface を**前層の語彙**に直した

| 旧 | 新 |
|---|---|
| `modTensor`(前層テンソル → 層化) | ★`preTensor`(前層テンソルのみ) |
| `IsInvertibleSheaf` | ★`IsInvertiblePre` |
| `sheafOf : Pic X → X.Modules` | ★`carrier : Pic X → X.PresheafOfModules` |

★★**局所自明な前層は自動的に層**なので「可逆前層の同型類」=「可逆層の同型類」であり、
**退化封じの強さは変わらない**。★★★これで `Found` 側と対応が付いた:

| Interface の要求 | `Found/Arakelov/` |
|---|---|
| `carrier` / `carrier_invertible` | `InvertiblePresheaf`(第 16) |
| `carrier_one` / `carrier_mul` | ★**rfl**(`one_carrier` / `mul_carrier`) |
| `carrier_injective` / `_surjective` | `PicPre` の商(第 17) |
| `group` | ★`CommGroup (PicPre X)`(第 17) |

### ★★★★残り 2 点の材料(2026-08-17 実測)

#### (1) `equivPicRing` —— **材料はほぼ揃っている**

`AlgebraicGeometry/Modules/Tilde.lean` に:

| 在庫 | 内容 |
|---|---|
| `tilde.adjunction.fullyFaithfulLOfIsIsoUnit` | ★★**`tilde` は充満忠実**(`IsIso unit` が instance) |
| `isIso_fromTildeΓ_iff` | `IsIso fromTildeΓ ↔ essImage` |
| `isIso_fromTildeΓ_of_presentation` | ★**表示を持てば `fromTildeΓ` は同型** |
| `instance : (tilde M).IsQuasicoherent` | `tilde` の像は準連接 |
| `LocallyFree.lean` の `IsLocallyFree → IsQuasicoherent` | ★**局所自由 ⟹ 準連接**(priority 100) |

★★★したがって「`Spec R` 上の可逆層は `tilde M` の形」は**出せる**。
★★残るのは **`tilde` がモノイダルであること**(`tilde (M ⊗_R N) ≅ tilde M ⊗ tilde N`)
——これは mathlib に無い。

#### (2) 引き戻し `f^* : Pic Y → Pic X`

★mathlib は `Scheme.Modules.pullback f` / `pullbackId` / `pullbackComp` /
`pullbackObjFreeIso`(`f^*(free I) ≅ free I`)を持つ。
★★★無いのは **`f^*` のモノイダル性**と**局所自明性の保存**。

### ★★★★★残り 2 点は同じ形の欠落である

    (1) tilde がモノイダル
    (2) f^* がモノイダル

★★どちらも「**関手がテンソル積を保つ**」——mathlib の
`Functor.Monoidal` / `CoreMonoidal` を作る仕事である。
★`pushforward₀OfCommRingCat` については mathlib が既に作っている
(`PushforwardZeroMonoidal.lean`)ので、**同じ型の仕事**であり見通しは立っている。


---

## §9-33 (2026-08-17) ★★★★★★B1 —— 20 ブロック。`equivPicRing` の乗法性が無料で出た

### ★★★★★★「前層の段で組む」判断の 2 つ目の配当

前層テンソルは**各開集合ごと**である(`tensorObj_obj` は定義そのもの):

    (F ⊗ G)(U) = F(U) ⊗_{𝒪(U)} G(U)

★★したがって `Γ = (·)(⊤)` は**乗法を保つ**。
★★★`picToSectionsPic X : PicPre X →* CommRing.Pic Γ(X, ⊤)` が取れた。

★★★★**前 turn に「`tilde` のモノイダル性が要る」と判定したのは誤りだった**
——見積り外し **6 度目**。

### ★機構(すべて mathlib の在庫)

| 何 | 道具 |
|---|---|
| 可逆加群であること | `L.isInv` を `⊤` で評価 → `Module.Invertible.left` |
| 乗法性 | `CommRing.Pic.mk_tensor` |
| well-defined | `CommRing.Pic.mk_eq_mk_iff` |
| `map_one` | `CommRing.Pic.mk_eq_one_iff` |

### ★★★★B1 の到達状況 —— `PicardData` の 10 フィールド

| フィールド | 状態 |
|---|---|
| `Pic` | ✅ `PicPre X`(第 17) |
| `group` | ✅ `CommGroup (PicPre X)`(第 17) |
| `carrier` / `carrier_invertible` | ✅ 第 16・19 |
| `carrier_one` / `carrier_mul` | ✅ **rfl**(第 16) |
| `carrier_injective` / `carrier_surjective` | ✅ 第 17・19 |
| `pullback` / `pullback_mul` / `pullback_id` / `pullback_comp` | ★残 |
| `equivPicRing` | ★★**乗法性は取得**(第 20)。残るは**全単射性** |

### ★★★残り 2 点の正体

#### (1) `equivPicRing` の全単射性 = **アフィンでの層と加群の対応**

★局所自明な前層 `F` on `Spec R` について `F ≅ tilde (F(⊤))` を言うこと。
★★材料は `Tilde.lean` に在る(`tilde` 充満忠実 /
`isIso_fromTildeΓ_of_presentation` / `IsLocallyFree → IsQuasicoherent`)が、
**我々の `Pic` は前層側**なので層側へ渡す段が要る。

#### (2) 引き戻し `f^*`

★mathlib の `Scheme.Modules.pullback` は**層**の側。
★★モノイダル性が未登録で、かつ前層側との橋が要る。

★★★**どちらも「前層 ↔ 層」の橋に帰着する。**
すなわち「**局所自明な前層は層である**」——これが B1 の最後の壁である。


---

## §9-34 (2026-08-17) ★★★★★★★**訂正: 「局所自明な前層は層」は偽だった**

### ★★★★★★何を間違えたか

第 16–20 ブロックで「**前層の段で `Pic` を組む**」道を採り、その根拠として

> 局所自明な前層は自動的に層である(層条件は局所的だから)

と繰り返し書いた。★★★**これは偽である。**

#### ★反例(ℙ¹ 上の直線束)

    X = ℙ¹,  F = 𝒪(1),  G = 𝒪(-1)

    (F ⊗_pre G)(ℙ¹) = H⁰(𝒪(1)) ⊗_k H⁰(𝒪(-1)) = k² ⊗_k 0 = 0
    𝒪(ℙ¹) = k

★前層テンソルは**各開集合ごと**(`tensorObj_obj`)なので、
`F ⊗_pre G` の大域切断は `0` であり `𝒪` と同型ではない。

#### ★★「層条件は局所的」も一般には偽

`F|_{U_i}` がすべて層でも `F` が層とは限らない
——`W` 上の貼り合わせを示すのに `{W ∩ U_i}` に沿う貼り合わせが要り、**循環する**。

### ★★★★何が壊れたか

| ブロック | 判定 |
|---|---|
| 1–15(層化テンソル・結合律・局所論法・局所自明性のテンソル閉性) | ✅ **無傷** |
| 16 `PicGroup`(`InvertiblePresheaf`) | ★★**`F ⊗_pre G ≅ 𝟙_` は直線束で成り立たない** |
| 17 `PicQuotient`(`PicPre` と `CommGroup`) | ★`PicPre` は Picard 群ではない |
| 18(Interface を前層語彙へ) | ★穴が再び開く(`carrier_surjective` が空回り) |
| 19–20(橋・`Γ` の乗法性) | ★誤った定義の上 |

★★★**4 ファイルを削除し、Interface を層の語彙(コミット `53c5868`)へ戻した。**

### ★★★★★教訓

★★★**「無い」と決める前に測る**を 6 度実践して成功したが、
今回は逆に**「有る(自動的に成り立つ)」と決めて測らなかった**。

★★★★**両方向に測る必要がある。**
「自動的に成り立つ」と思った補題こそ、**具体例で検算する**。
★今回は `ℙ¹` の `𝒪(1) ⊗ 𝒪(-1)` を 1 回計算すれば 5 ブロック前に気づけた。

### ★★正しい道(元の道)

`Pic X` は**層**の同型類でなければならない。したがって

    tensorModules L M := 層化 (L.val ⊗ M.val)     -- 第 2 ブロック

を乗法とし、第 13 ブロック(結合律)・第 14 ブロック(逆元)の上に組む。
★★★残る壁は変わらず **「層化が局所自明性を保つ」** ——
これは `(層化 P)|_V ≅ P|_V`(`P|_V` が層のとき)であり、
`Over V` 上の局所全単射から出す。★mathlib に該当補題は無い。


---

## §9-35 (2026-08-17) ★★★★★★壁の部品が**すべて**mathlib に在った

### ★訂正の直後に前進した

§9-34 で第 16–20 ブロックを撤回した(「局所自明な前層は層」が偽)。
★★正しい道(層側)に戻したうえで、残る壁

> **層化は局所自明性を保つ**:  `P|_V ≅ 𝟙_|_V`  ⟹  `(層化 P).val|_V ≅ 𝟙_|_V`

の**部品をすべて実測で確認した**。

### ★★★★★★4 部品すべて mathlib に在った

| 部品 | 在庫 | 状態 |
|---|---|---|
| `Over.forget V` が cocontinuous | `Sites/Over.lean` | ★取得(`isCocontinuous_overForget`) |
| 制限が局所全射を保つ | `Sites/PreservesLocallyBijective.lean` | ★取得(`isLocallySurjective_restrict`) |
| 制限が局所単射を保つ | 同上 | ★取得(`isLocallyInjective_restrict`) |
| **層どうしの局所全単射は同型** | `Sites/LocallyBijective.lean` の `Sheaf.isLocallyBijective_iff_isIso` | ★実測で確認 |

★★★★**「mathlib に該当補題なし」と 2 度判定したが、いずれも誤りだった**
——見積り外し **7 度目**。

### ★★★残る組み立て(機械的)

1. `η_P : P ⟶ (層化 P).val` は局所全単射(mathlib の instance、第 11 ブロックで確認済み)
2. 制限して `η_P|_V` も局所全単射(★本 turn で取得)
3. `P|_V` は層(`≅ 𝟙_|_V` かつ `𝟙_` の制限は層)
4. `(層化 P).val|_V` も層(制限は層を保つ——`op_comp_isSheaf_of_isSheaf`)
5. `Sheaf.isLocallyBijective_iff_isIso` で `η_P|_V` は同型
6. ⟹ `(層化 P).val|_V ≅ P|_V ≅ 𝟙_|_V`

★★**唯一の手間**は 3–5 で `PresheafOfModules` の射を `Sheaf J A` の射として
包み直すこと(`SheafOfModules.toSheaf` が `ReflectsIsomorphisms` を持つ)。

### ★★★★★この turn の 2 つの教訓

1. ★★★**「無い」と決める前に測る**(7 度成功)
2. ★★★★**「有る(自明に成り立つ)」と決める前にも測る**(1 度失敗)
   ——「局所自明な前層は層」を `ℙ¹` の `𝒪(1) ⊗ 𝒪(-1)` で検算していれば
   5 ブロック前に気づけた。**両方向に測ること。**


---

## §9-36 (2026-08-17) ★★★★★★★★B1 の中核が完成 —— `CommGroup (PicSheaf X)`

### ★取れたもの(`Found/Arakelov/` 15 ファイル、すべて sorry 0)

| 段 | 定理 |
|---|---|
| 壁 | ★`isLocallyTrivial_sheafify`(層化は局所自明性を保つ) |
| 閉性 | ★`isLocallyTrivial_tensorModules`(乗法が閉じる) |
| 並べ替え | ★`tensorRearrange`(`(A⊗B)⊗(A'⊗B') ≅ (A⊗A')⊗(B⊗B')`) |
| 型 | `InvSheaf X` / `one` / `symm` / `mul` |
| 商 | `PicSheaf X` |
| **群** | ★★★`instance : CommGroup (PicSheaf X)` |

### ★★★★★壁の破り方——4 部品すべて mathlib に在った

    η_P|_V : P|_V ⟶ (層化 P).val|_V
      ・局所全単射(`Sites/PreservesLocallyBijective.lean`)
      ・両側とも層(`Functor.op_comp_isSheaf` / `sheafCompose`)
    ⟹ 同型(`Sheaf.isLocallyBijective_iff_isIso`)
    ⟹ 前層加群の射としても同型(`toPresheaf` が同型を反映)

★★**「mathlib に該当補題なし」と 2 度判定したが、いずれも誤りだった。**

### ★★B1 の残り 2 点 —— どちらも「関手がテンソル積を保つ」

| 何 | 欠けているもの |
|---|---|
| 引き戻し `f^*` | `Scheme.Modules.pullback` は在るが**モノイダル性が未登録** |
| `equivPicRing` | `tilde` の**モノイダル性**(または `Γ` の準連接層上での乗法性) |

★mathlib は `pushforward₀OfCommRingCat` については既に作っている
(`PushforwardZeroMonoidal.lean`)ので、**同じ型の仕事**である。

### ★★★★★この session の方法論(記録)

| 姿勢 | 実績 |
|---|---|
| 「**無い**」と決める前に測る | ★**7 度成功** |
| 「**有る(自明)**」と決める前にも測る | ★★★**1 度失敗**(`ℙ¹` の `𝒪(1)⊗𝒪(-1)` で検算していれば防げた) |

★★★★**両方向に測る。**「自動的に成り立つ」と思った補題こそ具体例で検算する。


---

## §9-37 (2026-08-17) ★★★★★引き戻しのモノイダル性への道(実測)

### ★★mathlib の在庫

| 在庫 | 場所 |
|---|---|
| `(extendScalars f).Monoidal` | `ModuleCat/Monoidal/Adjunction.lean`(加群の係数変換) |
| `(restrictScalars f).LaxMonoidal` | 同上(随伴の右側から) |
| `Adjunction.rightAdjointLaxMonoidal` | ★`Monoidal/Functor.lean` 908 |
| `Adjunction.leftAdjointOplaxMonoidal` | ★★`Monoidal/Functor.lean` 1026 |
| `pushforward₀OfCommRingCat.Monoidal` | `Presheaf/PushforwardZeroMonoidal.lean`(**11 行**、μ・ε が恒等) |

### ★★★道筋

`PresheafOfModules.pullback φ := (pushforward φ).leftAdjoint` である。したがって

1. `(pushforward φ).LaxMonoidal` を示す
2. `Adjunction.leftAdjointOplaxMonoidal` で `pullback` が oplax monoidal
3. 比較射が同型であることを示して strong monoidal に上げる

★★これで `f^*(L ⊗ M) ≅ f^*L ⊗ f^*M` が出て、`PicardData.pullback` が書ける。

### ★同じ道が `equivPicRing` にも効く可能性

`tilde` も随伴の左側(`tilde ⊣ Γ`)なので、同じ手が使えるかもしれない。
★`Tilde.lean` の `tilde.adjunction` が材料。


---

## §9-38 (2026-08-17) ★★★引き戻しの残作業が 1 本に確定した

### ★実測(2026-08-17)

| 何 | 状態 |
|---|---|
| `ModuleCat.restrictScalars f` が lax monoidal | ★★**mathlib に在る**(`instLaxMonoidalRestrictScalars`) |
| その `μ` を取り出せる | ★探りで確認(`Functor.LaxMonoidal.μ`) |
| `PresheafOfModules.restrictScalars α` が lax monoidal | ★★★**無い** |
| `AlgebraicGeometry/Modules` にモノイダル性 | ★★★**無い** |

### ★★★残作業は 1 本

    (PresheafOfModules.restrictScalars α).LaxMonoidal

★`restrictScalars` は**成分ごと**の定義(`restrictScalarsObj` の `obj := fun X ↦ ...`)なので、
`ε` と `μ` は成分ごとに `ModuleCat` の在庫から取れる。
★★残るのは**自然性と 5 つの coherence** を前層レベルで書くこと。

★★★これが出れば:

    (pushforward φ).LaxMonoidal            -- pushforward₀(strong)⋙ restrictScalars(lax)
    → Adjunction.leftAdjointOplaxMonoidal  -- ★mathlib
    → pullback が oplax monoidal
    → 比較射が同型(可逆層なら局所的に恒等——第 15・16 ブロックが効く)
    → f^*(L ⊗ M) ≅ f^*L ⊗ f^*M
    → PicardData.pullback と pullback_mul

★★同じ随伴の手が `equivPicRing`(`tilde ⊣ Γ`)にも効く見込みである。

### ★規模

mathlib の `PushforwardZeroMonoidal.lean` は 11 行だが、あれは `μ`・`ε` が恒等の場合。
★`restrictScalars` の `μ` は**同型でない標準写像** `M ⊗_R N → M ⊗_S N` なので、
coherence を実際に書く必要がある。**mathlib PR 規模(100–200 行)**。


---

## §9-39 (2026-08-17) ★★★★★★★引き戻しの連鎖が繋がった

### ★取れたもの(第 18–20 ブロック、すべて sorry 0)

| # | 定理 | 備考 |
|---|---|---|
| 18 | `(resScalars α).LaxMonoidal` | ★★★**mathlib に無かった**。7 フィールド手書き |
| 19 | `(pushCR F R S α).LaxMonoidal` | ★strong ⋙ lax、**1 行** |
| 20 | `(PresheafOfModules.pullback (alphaR ..)).OplaxMonoidal` | ★`Adjunction.leftAdjointOplaxMonoidal` |

★★★**1 本(第 18)を作ったことで 3 段が自動で繋がった。**

### ★★第 18 ブロックの recipe(7 フィールドすべてに効いた)

1. `PresheafOfModules.hom_ext` で成分に落とす
2. `ModuleCat.hom_ext` + `TensorProduct.ext'` / `ext_threefold` で純テンソルへ
3. ★★**`show` で両辺を明示的に書き下す**(片側だけでは駄目)
4. `ModuleCat.restrictScalars_μ_tmul` / `restrictScalars_η` で書き換え
5. `rfl`

★`cat_disch` は 5 条すべてで落ちなかったが、この recipe は**5 条すべてで一発**だった。

### ★★残り —— 比較射 `δ` が同型であること

`OplaxMonoidal` は `δ : F(X ⊗ Y) ⟶ F X ⊗ F Y` を与える。
★★可逆層なら局所的に恒等なので同型になる——第 15・16 ブロック
(局所自明性のテンソル閉性・層化での保存)と第 8 ブロック(局所性)が効く。

★★★それが出れば `f^*(L ⊗ M) ≅ f^*L ⊗ f^*M` となり
`PicardData.pullback` / `pullback_mul` / `pullback_id` / `pullback_comp` が書ける
(後者 2 つは mathlib の `pullbackId` / `pullbackComp` から)。

### ★実装の罠(この turn)

★`(F.op ⋙ R) ⋙ forget₂` と `F.op ⋙ (R ⋙ forget₂)` は **defeq だが構文が違う**
——型注釈つきの別名を置くと暗黙引数が決まる。
★★`simp` は `Functor.LaxMonoidal.μ` を展開しない——`show` で適用形にする。


---

## §9-40 (2026-08-17) ★★★引き戻しの最後の 1 歩

### ★mathlib の在庫(実測)

    Functor.Monoidal.ofOplaxMonoidal [F.OplaxMonoidal] [IsIso (η F)] [∀ X Y, IsIso (δ F X Y)]

★★**oplax から strong へ上げる道具が在る**(`Monoidal/Functor.lean` 731)。

### ★★★ただし `δ` が全対象で同型なのは**偽**

★一般の加群では `f^*(M ⊗ N) → f^*M ⊗ f^*N` は同型でない(平坦性が要る)。
★★**可逆層に限れば**局所的に恒等なので同型になる。

### ★★残る 1 歩

    theorem isIso_delta_of_trivial (L M) (hL hM : IsLocallyTrivial ..) :
        IsIso (Functor.OplaxMonoidal.δ (pullback ..) L M)

★機構は第 8(局所性)・第 15(テンソル閉性)・第 16(層化での保存)ブロック
——「局所的に同型なら同型」を `δ` に当てる。

★★★これが出れば `f^*(L ⊗ M) ≅ f^*L ⊗ f^*M` となり、
`PicardData` の `pullback` / `pullback_mul` が書ける
(`pullback_id` / `pullback_comp` は mathlib の `pullbackId` / `pullbackComp` から)。


---

## §9-41 (2026-08-18) ★★★★★★★★§9-40 の訂正 —— `δ` は**全対象で同型**である

### ★★★誤っていた記述(§9-40)

> ★★★ただし `δ` が全対象で同型なのは**偽**
> ★一般の加群では `f^*(M ⊗ N) → f^*M ⊗ f^*N` は同型でない(平坦性が要る)。

★★★★**これは誤りである。**平坦性が要るのは `f^*` の**完全性**であって、
**モノイダル性ではない**。

### ★★★★正しい事実

**`f^*` は 𝒪 加群について strong monoidal である**(Stacks 01BJ / Hartshorne II.5)。

    f^*F = 𝒪_X ⊗_{f⁻¹𝒪_Y} f⁻¹F

    f^*(F ⊗_{𝒪_Y} G)
      = 𝒪_X ⊗_{f⁻¹𝒪_Y} f⁻¹(F ⊗_{𝒪_Y} G)
      = 𝒪_X ⊗_{f⁻¹𝒪_Y} (f⁻¹F ⊗_{f⁻¹𝒪_Y} f⁻¹G)   ★f⁻¹ が strong monoidal
      = (𝒪_X ⊗ f⁻¹F) ⊗_{𝒪_X} (𝒪_X ⊗ f⁻¹G)         ★係数拡大が strong monoidal
      = f^*F ⊗_{𝒪_X} f^*G

★★2 段とも mathlib に対応物がある:
`(extendScalars f).Monoidal` は**在る**(`ModuleCat/ChangeOfRings.lean`)。

### ★★★★★これが意味すること

★可逆層に限る必要が**無い**——`δ` は**恒等的に同型**である。
★★したがって `PicardData.pullback` の障害は「可逆性の局所論法」ではなく、
**「Lean で strong monoidal 性をどう出すか」**という技術問題に変わった。

### ★★Lean での道(3 案、いずれも道であって壁ではない)

| 案 | 中身 | 難度 |
|---|---|---|
| (a) 生成元 | `δ` は両辺とも各変数で余極限を保ち、`free (yoneda X)` 上で一致 | ★★中(mathlib に `isColimitFreeYonedaCoproductsCokernelCofork` が在る) |
| (b) 局所論法 | 第 8・15・16 ブロックで「局所同型 ⟹ 層化後に同型」 | ★★★`f^*` の局所計算が要る |
| (c) mate | `Adjunction.IsMonoidal` を作る | ★★★★coherence 4 条 |

★★★**(a) が本命**である。`free` 上では mathlib の
`SheafOfModules.pullbackObjFreeIso` が既に同型を与えている。

### ★★★★★★方法論の再確認

★★§9-40 で私は「平坦性が要る」と**測らずに書いた**。
★★★これは第 16–20 ブロックで一度踏んだ失敗
(「局所自明な前層は自動的に層」——偽)と**同じ型**である。

★**「無い」と決める前に測る**(7 回成功)
★★**「有る」と決める前にも測る**(1 回失敗)
★★★**「偽」と決める前にも測る**(本件、2 回目の失敗)

——3 つとも同じ規則の系である: **決める前に測る。**


---

## §9-42 (2026-08-18) ★★★★★★第 21–24 ブロックと、`δ` 同型化の残り 1 点

### ★★★★★★この turn で取れたもの(すべて sorry 0、push 済み)

| # | 定理 | 一言 |
|---|---|---|
| 21 | `opensMapFinal` / `pullbackUnitIso` | ★★★★**構造層の引き戻しは構造層** |
| 22 | `pullbackDelta` / `sheafifyPullbackApp` | ★★★★**比較射 δ がスキームの射に効く**(`alphaR` と `rfl` で一致) |
| 23 | 余極限保存 3 本 / `opensMap_inf` | ★生成元による道の前提 |
| 24 | `pullbackFreeYonedaIso` | ★★★★★**抽象な引き戻しの具体値**(余表現の一意性) |

★★★**`pullback_id` / `pullback_comp` は mathlib に在る**
(`Scheme.Modules.pullbackId` / `pullbackComp`)。`PicSheaf` は同型による商なので、
自然同型がそのまま等式になる。

### ★★残り 1 点 —— `free (yoneda V) ⊗ free (yoneda W) ≅ free (yoneda (V ⊓ W))`

★構成そのものは**書けた**(`app` は型検査を通る):

    app Z := ModuleCat.freeDesc (↾fun q => freeMk q.1 ⊗ₜ freeMk q.2)

★★**通らないのは naturality の証明**である。

### ★★★★★★障害の正体(2026-08-18 実測)——**`Ring` インスタンスの二重路**

`↑((R ⋙ forget₂ CommRingCat RingCat).obj Z)` には `Ring` が **2 経路**で付く:

| 経路 | 出どころ |
|---|---|
| (a) | `RingCat` の構造(`freeObj` / `restrictScalars` が使う) |
| (b) | `CommRing.toRing`(`Presheaf/Monoidal.lean` 34 行の `CommRing` インスタンス経由。テンソルが使う) |

★★**(a) と (b) は defeq だが構文が違う。**
そのため `simp` / `rw` は `ConcreteCategory.hom (f ≫ g) y` **にすら当たらない**
——`CategoryTheory.comp_apply` は `@[simp]` なのに発火しない。

★★★これは第 11 ブロックで記録した
「インスタンス束縛子は型の書き方の違いをまたげない」の**兄弟**である。
違いは、あちらが**探索の失敗**で、こちらが**書き換えの失敗**であること。

### ★塞ぎ方の候補(次の turn で測る)

1. `letI` で `Ring` インスタンスを片方に固定してから構成する
2. `ModuleCat.of` で束ね直し、両辺を同じ経路の項にする
3. naturality を `LinearMap` の段まで降ろし、`TensorProduct.ext'` で純テンソルにしてから
   `show` で両辺を書き下す(第 18 ブロックの recipe と同じ手口)

★★**3 が本命**である——第 18 ブロックでは同種の壁を
「`show` で両辺を明示的に書き下す」で 7 フィールドすべて抜けた。


---

## §9-43 (2026-08-18) ★★★★★★第 25・26 ブロックと、残り 2 点の見取り図

### ★★★★★★二重路を抜けた(第 25 ブロック)

`Ring` インスタンスの二重路((a) `RingCat` 由来 / (b) `CommRing.toRing` 由来)は
**2 つの組み合わせ**で抜けた:

1. **加群の段で naturality を書かない**——`freeObjDesc` に渡せば mathlib が担う。
   ★型前層の段なら等式は**要素の等式**である。
2. **`erw` を使う**——`rw` は `instances` 透明度で照合するので落ちるが、
   `erw` は `default` 透明度なので**二重路を越えられる**。

★★★これは記録に値する: 「`simp` が `@[simp]` 付きの補題にすら当たらない」ときは
**インスタンス経路を疑い、`erw` を試す**。

### ★★生成元の段が 4 本で閉じた

    第 23: 逆像は ⊓ を保つ                        `opensMap_inf`
    第 24: f^*(free (yoneda V)) ≅ free (yoneda (f⁻¹V))
    第 25: free (F ⊗ G) ≅ free F ⊗ free G
    第 26: yoneda V ⊗ yoneda W ≅ yoneda (V ⊓ W)

★★★**両辺の対象が生成元の上で同型である**ことが確定した:

    f^*(free(yoneda V) ⊗ free(yoneda W)) ≅ free(yoneda (f⁻¹V ⊓ f⁻¹W))
                                         ≅ f^*free(yoneda V) ⊗ f^*free(yoneda W)

### ★★★★残り 2 点

| # | 主張 | 材料 |
|---|---|---|
| (a) | `δ` **自身**が上の同型と一致すること | ★mathlib は `δ X Y = adj.homEquiv.symm ((unit ⊗ₘ unit) ≫ μ)` と**具体的に**定義している。★★`pushforwardCompCoyonedaFreeYonedaCorepresentableBy` の `homEquiv` も `freeYonedaEquiv` で**具体的**である |
| (b) | 余極限による持ち上げ | ★第 23 の余極限保存 3 本 + `isColimitFreeYonedaCoproductsCokernelCofork` |

### ★★★★★★重要な観察 —— `δ` を経由しなくてもよい可能性

`PicSheaf` は**同型による商**である。したがって `pullback_mul` に要るのは
**対象の同型** `f^*(L ⊗ M) ≅ f^*L ⊗ f^*M` だけであり、
★★**それが `δ` であることは要らない**。

★★★これは (a) を迂回できる可能性を示す——ただし対象の同型を出すにも
結局は生成元と余極限を通す必要があり、(b) は残る。


---

## §9-44 (2026-08-18) ★★★★★★同じ型の穴が 4 例——一般規則を抽出した

### ★★4 例

| # | 場所 | 穴 | 状態 |
|---|---|---|---|
| 1 | C1 `evalAffine` | 値がアフィンでしか固定されていなかった | ★塞いだ(2026-08-17) |
| 2 | B1 `Pic` | 非アフィンで自由だった | ★塞いだ(`sheafOf` 追加、2026-08-17) |
| 3 | B1 `pullback` | `sheafOf` と結ばれていなかった | ★塞いだ(`sheafOf_pullback` 追加、2026-08-18) |
| 4 | G7 `omega` | **曲線と結ばれていない** | ★★**未**(検査だけ書いた) |

★どれも `Check/` で**機械的に**穴の実在を示した。

### ★★★★★★一般規則

> **`→ Type` の posit(台を仮定する欄)は、必ず
> 「その欄と入力データの両方を言及する条件」を 1 本以上持たねばならない。**

★★階数・濃度・群構造だけでは足りない——それらは**入力を無視する定数**で満たせる。

| 足りない条件の例 | 通ってしまう witness |
|---|---|
| 「階数 1」だけ | `omega L W := 𝓞 L`(曲線を無視) |
| 「群である」だけ | `Pic X := Multiplicative ℤ`(スキームを無視) |
| 「アフィンで一致」だけ | 非アフィンで自由 |

### ★★★監査すべき残り(`Interface/` の `→ Type` posit、2026-08-18 実測で 20 件超)

    APic.APic / APic.AlgPoint / ArcSpace.Metric
    Torsion.tateModule
    GenEll/AbcSetup.Curve / GaloisRep.EllClass / GaloisRep.Curve
    HeightTheory.Point / ABundle / Divisor / GenericClass
    TateLocal.LocalField / Curve / StableLine / MultExt
    IUTchIII/PilotObjects.Amb
    NCBelyi/BelyiSetup.Curve / Point / ProjPoint / Mor / NumField

★★★**「達成」を数える前に、この監査を通すこと。**
★通っていない欄が残ったまま埋めると、**誤った達成**が数えられる。


---

## §9-45 (2026-08-18) ★★★★★★第 25–33 ブロック——`δ` の同型性が「生成元の 1 点」に絞られた

### ★★★★★★連鎖の全体像

| # | 取れたもの | 一言 |
|---|---|---|
| 25 | `free (F ⊗ G) ≅ free F ⊗ free G` | ★`Ring` インスタンスの二重路を `erw` で抜けた |
| 26 | `yoneda V ⊗ yoneda W ≅ yoneda (V ⊓ W)` | ★前順序だから表現可能 |
| 27 | `f^*(free(yV) ⊗ free(yW)) ≅ f^*free(yV) ⊗ f^*free(yW)` | ★★**対象の同型** |
| 28 | `isIso_app_of_isColimit` | ★★余極限で同型を持ち上げる(**mathlib に無し**) |
| 29 | `isIso_app_of_freeYoneda` | ★★生成元で同型なら全対象で同型 |
| 30 | `isIso_pullbackDelta_of_free` | ★★★★**二重の持ち上げ** |
| 31 | `pullbackUnit_app_free` | ★★★抽象な `unit` の書き下し |
| 32 | `pushforward_μ_app_tmul` | ★★`μ` は切断の上で恒等 |
| 33 | `homEquiv_pullbackDelta` / `mathlibCorep_homEquiv` | ★計算の手がかり |

### ★★★★残り 1 点

    IsIso (δ (free (yoneda V)) (free (yoneda W)))

★これだけを示せば、第 30 ブロックの二重持ち上げで**全対象**に伝播する。
★★材料は揃った:

    δ の adjunct = (unit ⊗ₘ unit) ≫ μ        (第 33)
    unit は freeYonedaEquiv で書ける          (第 31 + 第 33)
    μ は純テンソルを純テンソルに送る            (第 32)

### ★★★実装の罠(この turn で 3 つ)

1. **`Ring` インスタンスの二重路** → `erw`(`default` 透明度)で抜ける
2. **合成の余極限保存が探索できない** → 汎用の圏で補題を書く
3. **保存性をインスタンス束縛子で受けると渡せない** → `(hA hB : ..)` と明示引数にする

★★2 と 3 は第 11 ブロックの規則(型の書き方をまたげない)の 6・7 例目である。


---

## §9-46 (2026-08-18) ★★★★★★★最終段の骨組みが型検査を通った——残りは要素の等式 1 本

### ★★★★★★示すべき等式(型検査済み、2026-08-18)

    (pullbackPre f).map (iotaY V W) ≫ pullbackDelta f (freeY V) (freeY W)
      = (pullbackFreeYonedaIso f (V ⊓ W)).hom
        ≫ (targetIso f V W).inv
        ≫ ((pullbackFreeYonedaIso f V).symm ⊗ᵢ (pullbackFreeYonedaIso f W).symm).hom

ここで

    iotaY V W  := (freeYonedaInfIso (R := Y.presheaf) V W).inv        -- 第 26
    targetIso  := freeYonedaInfIso (R := X.presheaf) (f⁻¹V) (f⁻¹W)     -- 第 26(X 側)

★★★**右辺は同型の合成である。**したがってこの等式が証明できれば

    IsIso ((pullbackPre f).map (iotaY V W) ≫ δ)

が出て、`(pullbackPre f).map (iotaY V W)` が同型だから **`IsIso δ`** が出る。
★第 30 ブロックの二重持ち上げで**全対象**に伝播する。

### ★★証明の骨組み(ここまで通っている)

    refine (adj.homEquiv _ _).injective ?_          -- ★随伴の単射性
    erw [Adjunction.homEquiv_naturality_left]       -- ★★`erw` が要る(型の書き方の違い)
    rw [homEquiv_pullbackDelta]                     -- ★第 33: δ の adjunct = (unit ⊗ₘ unit) ≫ μ
    refine PresheafOfModules.freeYonedaEquiv.injective ?_   -- ★生成元だから元で決まる
    simp only [PresheafOfModules.freeYonedaEquiv_comp]      -- ★合成則
    -- ★★★残り: 要素の等式 1 本

### ★★★残った要素の等式の中身(数学としては確認済み)

左辺の元は

    unit_P(freeMk i_V) ⊗ₜ unit_Q(freeMk i_W)        (第 31 + 第 32 で計算できる)

であり、右辺は `freeMk (𝟙 U₀)` の像である。★`U₀ = f⁻¹V ⊓ f⁻¹W = f⁻¹(V ⊓ W)`(第 23、`rfl`)。

★★**両者が一致することは、`unit` が生成元を生成元に送ることに他ならない。**

### ★注意

★この骨組みは `sorry` を含むので `Found/` には置けない(規則どおり)。
★★本節に全文を記録した——次の turn はここから再開する。


---

## §9-47 (2026-08-18) ★★★★★★★★最終段の証明が数学として完成した

### ★★★★★★示すべき等式(第 30 ブロックで全対象に伝播する)

    f^*(iotaY V W) ≫ δ = (pullbackFreeYonedaIso f (V ⊓ W)).hom
                          ≫ (targetIso f V W).inv
                          ≫ ((pullbackFreeYonedaIso f V).symm ⊗ᵢ (..W).symm).hom

★右辺は同型の合成なので、これで `IsIso δ` が出る。

### ★★両辺を要素まで降ろすと

**左辺**(`adj.homEquiv` → `freeYonedaEquiv` → 合成則):

    unit_P.app (op (V ⊓ W)) (freeMk i_V) ⊗ₜ unit_Q.app (op (V ⊓ W)) (freeMk i_W)

- `freeYonedaEquiv (iotaY V W) = freeMk i_V ⊗ₜ freeMk i_W`
  (第 26 の `yonedaInfIso.inv (𝟙) = (i_V, i_W)` + 第 25 の `freeTensorHom_app_freeMk`)
- `(unit_P ⊗ₘ unit_Q)` は純テンソルを成分ごとに送る
- `μ` は純テンソルを純テンソルに送る(第 32)

**右辺**(`homEquiv_naturality_right` → 第 34 → 第 35):

    inv_V.app (op U₀) (freeMk j_V) ⊗ₜ inv_W.app (op U₀) (freeMk j_W)

- `freeYonedaEquiv (adj.homEquiv hom) = hom.app (freeYonedaEquiv unit) = freeMk (𝟙 U₀)`
  ——★★これが第 35 ブロック(`unit` は生成元を生成元に送る)
- `targetIso.inv (freeMk 𝟙) = freeMk j_V ⊗ₜ freeMk j_W`(X 側で同じ計算)
- `(inv_V ⊗ₘ inv_W)` を当てる

### ★★★★★★一致する理由(1 行)

> **両辺とも `freeYonedaEquiv (unit)` の `U₀` への制限である。**

★左辺: `freeMk i_V = P.map i_V.op (freeMk 𝟙)` なので、`unit_P` の自然性から
`unit_P.app (freeMk i_V)` は `unit_P.app (freeMk 𝟙)` の制限。
★★右辺: `inv_V` の自然性から `inv_V.app (freeMk j_V)` も同じ元の制限。
★★★そして第 34 ブロックが `freeYonedaEquiv (unit) = freeYonedaEquiv (inv)` を与える。

### ★残る作業

★Lean での鎖(`NatTrans.naturality_apply` を 2 回、その他は既存の補題)。
★★**数学的な不確定性は無い**——記号の運搬だけである。


---

## §9-48 (2026-08-18) ★★★★★★★★★★`δ` が同型になった —— 引き戻しは strong monoidal

### ★★★★★★★★取れたもの(第 40・41 ブロック、sorry 0)

    isIso_pullbackDelta   : ∀ P Q, IsIso (δ : f^*(P ⊗ Q) ⟶ f^*P ⊗ f^*Q)      ★前層の段
    pullbackSheafTensorIso: f^*(L ⊗ M) ≅ f^*L ⊗ f^*M                        ★層の段

★★★これは **mathlib に無い**(2026-08-18 実測)。平坦性を要しない事実(Stacks 01BJ)。

### ★★★★道のり —— 23 ブロック(第 18–40)

| 段 | ブロック | 内容 |
|---|---|---|
| モノイダル | 18–22 | 係数変換 lax → 押し出し lax → 引き戻し oplax → **比較射 δ** |
| 生成元 | 23–27 | 余極限保存・生成元の具体形・`⊓` の保存 |
| 持ち上げ | 28–30 | 余極限で同型を持ち上げ、δ を生成元 1 点に帰着 |
| 計算 | 31–39 | **抽象な `unit` を生成元の上で計算し切る** |
| 完成 | 40–41 | δ が同型 → 層の段へ |

### ★★★★★★一番難しかったところ

`PresheafOfModules.pullback` は**随伴関手定理**で作られており、
`unit` に具体的な記述が**無い**。

★★★抜け方は `CorepresentableBy` の一意性(第 24 ブロック)である:

    同じ関手を余表現する対象は同型
      ⟹ f^*(free (yoneda V)) ≅ free (yoneda f⁻¹V)
      ⟹ その同型の `inv` が `unit` を表す(第 31)
      ⟹ `freeYonedaEquiv` で**要素**になる(第 33–35)

★★核心は 1 行: **両辺とも `freeYonedaEquiv (unit)` の制限である**(第 36・38)。

### ★★残り(B1)

| # | 残作業 | 見通し |
|---|---|---|
| (a) | `f^*` が局所自明性を保つ | ★`pullbackComp` を 2 回——開埋め込みなら底変換は要らない |
| (b) | `PicSheaf` の `pullback` と 3 公理 | ★★`pullback_id`/`comp` は mathlib から。`pullback_mul` は第 41 から |
| (c) | `equivPicRing`(アフィン比較) | ★★★`Tilde.lean` の `tilde` が充満忠実。未着手 |


---

## §9-49 (2026-08-18) ★★★★第 41 ブロックの仮定は**放置できない**——次に測る点

### ★★仮定の中身

第 41 ブロック `pullbackSheafTensorIso` は

    hM  : IsLocallyRankOne X ((pullbackPre f).obj M.val)
    hLs : IsLocallyRankOne X ((sheafifyFunctor X).obj ((pullbackPre f).obj L.val)).val

を要求する。★これらは第 11 ブロック(`sheafifyTensorRight` / `Left`)から来る。

### ★★★★問題 —— **前層の引き戻しの切断は階数 1 とは限らない**

`PresheafOfModules.pullback` は随伴関手定理で作られており、
切断は余極限である。★★**層化して初めて可逆層になる**のであって、
前層のままで局所階数 1 とは限らない。

★★★したがって `hM` は**そのままでは示せない可能性が高い**。

### ★★★測るべきこと(次の turn)

| # | 問い |
|---|---|
| 1 | `IsLocallyRankOne X (f^*_pre 𝒪_Y)` は成り立つか(第 21 ブロックの `η` は層の段の話) |
| 2 | 第 11 ブロックの仮定を**層の段**(`IsLocallyRankOne` の sheafify 版)に弱められるか |
| 3 | 第 9・10 ブロック(局所基底からの局所単射性)は前層の基底を本当に要するか |

★★★★**3 が本丸**である。局所全単射性は局所的な性質なので、
「層化すれば階数 1」で足りる可能性がある——ただし
テンソルは切断ごとなので、**測らないと分からない**。

### ★方法論

★★これは「**決める前に測る**」の適用である——
仮定を置いたまま先へ進むと、後で**討ち死に**する。
★§9-40 で「δ は全対象で同型でない」と測らずに書いて誤った
(§9-41 で訂正)のと同じ型の危険である。
