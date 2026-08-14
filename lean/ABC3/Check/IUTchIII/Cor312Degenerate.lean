import ABC3.Skeleton.IUTchIII.Cor312

/-!
# [IUTchIII] Corollary 3.12 — 不定性を自明化すると statement は自明になる

`Skeleton/IUTchIII/Cor312.lean` の `cor_3_12` について、
**明らかな退化**を機械にかける。

退化のさせ方は原文の非対称性が指す一点だけ:
Θ 側が q 側と違うのは「**可能な像**の和集合の**正則包**」を測る点なので、

- 可能な像を **`{qImage}` ただ1つ**にする(= 不定性が何も動かさない)
- 正則包を **恒等**にする

とすれば Θ 側と q 側は同じものを測ることになる。

## 結果(下で証明、`sorry` 無し)

`thetaLogVol = qLogVol`(`trivialised_thetaLogVol_eq_qLogVol`)。
したがって Corollary 3.12 の結論は**等号として成り立ち**、
しかも台 `A`・対数体積 `lv`・像 `S`・`|log(q)|` の値が**何であっても**成り立つ
(`trivialised_satisfies_cor_3_12`)。不等式は何も言っていない。

これは Scholze–Stix の指摘

> the critical [IUTT-3, Theorem 3.11] does not become false, but **trivial**

を、**我々のモデルの上で**再現したものである。

## ★これは原典への判定ではない(必読)

自明化したのは我々が置いた `Interface`(`PilotObjectData`)であって、原典ではない。
自明化が起きた原因は次のどれでもありうる:

1. **我々の `Interface` が弱い** — 原文が (Ind1)(Ind2)(Ind3) や正則包に課している
   条件を、我々が逐語から拾い切れていない。`PilotObjectData` は
   `hull` に何の条件も課しておらず、`possibleThetaImages` にも課していない。
2. **必要な数学が未構築** — 条件を書こうにも mono-theta 環境・log-shell・
   Frobenioid・tempered 基本群が無い(mathlib に 0 件)。
3. 原典側の飛躍。

**既定は 1 である**(PLAN §5-2)。3 を名乗るには複数の独立な型設計で同じ壁に当たることと
falsifier が要る。現時点でそれは満たしていないので、`Gap/` には登録しない。

**falsifier(何が起きればこの退化が 1 に落ちるか)**: 原文の逐語から
「`hull` は包含を保つ」「`possibleThetaImages` は空でない」等の条件を
**追加で読み取れた**とき、この退化 witness はそれを満たさなくなりうる。
そのときは本ファイルの witness を作り直す。

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.IUTchIII

open ABC3.Interface.IUTchIII ABC3.Skeleton.IUTchIII

variable {A : Type} (lv : Set A → WithTop ℝ) (S : Set A) (a : ℝ)

/-- **不定性を自明化した `PilotObjectData`**。

`Interface` が課している条件(`qAbs_pos` と `qLogVol_eq`)は残したまま、
Θ 側だけを潰す——可能な像は `{S}` の1つ、正則包は恒等。
台 `A`・対数体積 `lv`・像 `S`・`a = |log(q)|` は**任意**でよい。 -/
noncomputable def trivialised (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    PilotObjectData where
  Amb := A
  topology := τ
  logVol := lv
  hull := id
  possibleThetaImages := {S}
  qImage := S
  qAbs := a
  qAbs_pos := ha
  qLogVol_eq := h
  outputLogVolumes := {x : ℝ | (x : WithTop ℝ) ≤ lv S}
  outputLogVolumes_eq := by simp
  qLogVol_mem := le_of_eq h.symm
  thetaUnion_isCompact := by simpa using hS
  logVol_hull_ne_top_of_isCompact := hfin

/-- ★**退化の核心**: 自明化すると Θ 側と q 側は**同じ値**になる。 -/
theorem trivialised_thetaLogVol_eq_qLogVol (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    thetaLogVol (trivialised lv S a τ ha h hS hfin) = qLogVol (trivialised lv S a τ ha h hS hfin) := by
  simp [thetaLogVol, qLogVol, trivialised]

/-- ★**Corollary 3.12 の結論が、残りのデータが何であっても成り立つ**。

`sorry` 無し。すなわち自明化した `Interface` の下では、
`cor_3_12` は**内容を持たない**。 -/
theorem trivialised_satisfies_cor_3_12 (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) :
    thetaLogVol (trivialised lv S a τ ha h hS hfin) ≠ ⊤ ∧
      qLogVol (trivialised lv S a τ ha h hS hfin) ≤ thetaLogVol (trivialised lv S a τ ha h hS hfin) := by
  refine ⟨?_, le_of_eq (trivialised_thetaLogVol_eq_qLogVol lv S a τ ha h hS hfin).symm⟩
  rw [trivialised_thetaLogVol_eq_qLogVol lv S a τ ha h hS hfin, qLogVol]
  show lv S ≠ ⊤
  rw [h]
  exact WithTop.coe_ne_top

/-- 「i.e.」の形も自明に成り立つ——`C_Θ ≥ −1` は等号の場合 `C_Θ = −1` で達成される。
すなわち退化した witness は、**abc へ効く形の主張も**空虚に満たす。 -/
theorem trivialised_satisfies_CTheta (τ : TopologicalSpace A) (ha : 0 < a)
    (h : lv S = ((-a : ℝ) : WithTop ℝ)) (hS : @IsCompact A τ S)
    (hfin : ∀ U : Set A, @IsCompact A τ U → lv U ≠ ⊤) (C : ℝ) (hC : thetaLogVol (trivialised lv S a τ ha h hS hfin) ≤ ((C * a : ℝ) : WithTop ℝ)) :
    -1 ≤ C := by
  rw [trivialised_thetaLogVol_eq_qLogVol lv S a τ ha h hS hfin, qLogVol] at hC
  show (-1 : ℝ) ≤ C
  have h4 : -a ≤ C * a := by
    have : ((-a : ℝ) : WithTop ℝ) ≤ ((C * a : ℝ) : WithTop ℝ) := by
      rw [← h]; exact hC
    exact_mod_cast this
  nlinarith [ha]

-- 退化 witness であるにもかかわらず、受理ゲートを全部通ることの確認
#print axioms trivialised_satisfies_cor_3_12

/-! ## ★★ 反証は死んだ — 原文 p.184 と p.175 を足した後の最終測定(2026-08-14)

以前ここには **`cor_3_12` は現在の `Interface` の下では偽である**という定理
(`cor_3_12_refutable_under_current_interface`)と、その反例
(`possibleThetaImages := {∅}`、空集合の対数体積を `⊤` にしたもの)が置いてあった。

原文へ当たり直して**2段階**で足した結果、反例は**構成できなくなった**。

| 追加したもの | 出所 | 効いた先 |
|---|---|---|
| `outputLogVolumes` / `outputLogVolumes_eq` / `qLogVol_mem` | 物理 p.184、Step (xi-e)(xi-f) | **(b) 不等式** |
| `thetaUnion_isCompact` / `logVol_hull_ne_top_of_isCompact` | 物理 p.175 + p.31 + p.127 | **(a) 有限性** |

結果:

| | 内容 | 結果 |
|---|---|---|
| **(b)** | `qLogVol D ≤ thetaLogVol D` | **証明できる** |
| **(a)** | `thetaLogVol D ≠ ⊤` | **証明できる** |

すなわち `Skeleton.IUTchIII.cor_3_12` は **`sorry` 無しで通る**。

### ★★判定(課題の第3段)— **祝う話ではない**

**(a)(b) の両方が通ったので、我々のモデルの中で Corollary 3.12 は自明になった。**
これは Scholze–Stix の「trivial, not false」を我々のモデル上で再現したものである。

**我々は何も導出していない。** `cor_3_12` の証明は、`Interface` の
`logVol_hull_ne_top_of_isCompact ∘ thetaUnion_isCompact` と `qLogVol_mem` を
そのまま組み合わせただけである。(xi-e)(xi-f) も p.175 のコンパクト性も、
すべて**原文の証明から仮説として輸入**した。

**これは原典への判定ではない。** 原文は実際に証明を書いており、
我々が写していないのはその中身の方である。

### 旧反例はどの条件で死んだか(下で機械的に確定させる)

★`⋃₀ {∅} = ∅` は**コンパクト**なので、`thetaUnion_isCompact` は**満たしてしまう**。
死因は **`logVol_hull_ne_top_of_isCompact`**——p.31 の
「compact, hence of finite log-volume」に当たる橋の方である。

### ★残った問い: 原文は "easily verified" の検証を書いていない

`thetaUnion_isCompact` の原文の根拠は p.175 の「the **[easily verified]** compactness」
だけで、**検証は書かれていない**(2026-08-14 実測: 論文全体で `compact` は 33 件、
`¹,°𝒰_{j,v_ℚ}` のコンパクト性を確立する箇所は **0 件**)。方針としては
p.128 の「**all the indeterminacies that occur in the theory are compact** is in some
sense one important theme in the present series of papers」がある。
また `bounded family`(p.130、Remark 3.9.5, (vi))は**用語の定義**として1件あるのみで、
Theorem 3.11 の出力が bounded family であるとは述べていない。

したがって次に転写すべきは **Theorem 3.11, (i), (a) のモノ解析的整構造
`I(…) ⊆ I^ℚ(…)`(= log-shell、p.31 で compact と明示)と Proposition 3.9** である。
そこまで写して初めて `thetaUnion_isCompact` を**導出**できる可能性が出る。
-/

open Classical in
/-- 旧反例が使っていた対数体積:「空集合は `+∞`、それ以外は `−1`」。 -/
noncomputable def lvUnit (s : Set Unit) : WithTop ℝ :=
  if () ∈ s then ((-1 : ℝ) : WithTop ℝ) else ⊤

/-- ★**旧反例の死因を機械的に確定させる**。

`∅` はコンパクトなので `thetaUnion_isCompact` は満たされてしまう。
破れるのは `logVol_hull_ne_top_of_isCompact`(`hull = id` なので `logVol ∅ = ⊤`)。
すなわち反証を殺したのは **p.31 の「compact, hence of finite log-volume」の橋**である。 -/
theorem oldCounterexample_dies_on_the_bridge :
    IsCompact (∅ : Set Unit) ∧ ¬ (∀ U : Set Unit, IsCompact U → lvUnit (id U) ≠ ⊤) := by
  refine ⟨isCompact_empty, fun h => ?_⟩
  exact h ∅ isCompact_empty (by simp [lvUnit])

#print axioms oldCounterexample_dies_on_the_bridge

/-!
## 読み

`cor_3_12` の結論は、自明化した `PilotObjectData` の下で **等号**になる。
不等式 `≥` が主張になるのは、Θ 側の「可能な像」が本当に**複数**あり、
その和集合の正則包が q 側の像より**真に大きい**ときだけである。
その「複数性」を作っているのが (Ind1)(Ind2)(Ind3) だから、
**不定性の内容を書き下さない限り、Corollary 3.12 は内容を持たない。**

これは `Check/PGC/InertiaDegeneracy.lean`(惰性群を自由なデータに置くと
`I_K := ⊥` でも `⊤` でも Corollary 1.3 が成り立つ)と同型の失敗であり、
`Meta/Calibration.lean` が実演した「型は付くが何も言っていない」の3例目である。

**次に何をすべきか**: `possibleThetaImages` と `hull` に、
原文の逐語から**追加で読み取れる条件**があるかを洗う。
無ければ、それは「原文がこの段で条件を与えていない」ことの実測になり、
`Gap/` の候補として §5-2 のトリアージにかける(既定は①、我々のモデル化の誤り)。
-/

end ABC3.Check.IUTchIII
