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
noncomputable def trivialised (ha : 0 < a) (h : lv S = ((-a : ℝ) : WithTop ℝ)) :
    PilotObjectData where
  Amb := A
  logVol := lv
  hull := id
  possibleThetaImages := {S}
  qImage := S
  qAbs := a
  qAbs_pos := ha
  qLogVol_eq := h

/-- ★**退化の核心**: 自明化すると Θ 側と q 側は**同じ値**になる。 -/
theorem trivialised_thetaLogVol_eq_qLogVol (ha : 0 < a) (h : lv S = ((-a : ℝ) : WithTop ℝ)) :
    thetaLogVol (trivialised lv S a ha h) = qLogVol (trivialised lv S a ha h) := by
  simp [thetaLogVol, qLogVol, trivialised]

/-- ★**Corollary 3.12 の結論が、残りのデータが何であっても成り立つ**。

`sorry` 無し。すなわち自明化した `Interface` の下では、
`cor_3_12` は**内容を持たない**。 -/
theorem trivialised_satisfies_cor_3_12 (ha : 0 < a) (h : lv S = ((-a : ℝ) : WithTop ℝ)) :
    thetaLogVol (trivialised lv S a ha h) ≠ ⊤ ∧
      qLogVol (trivialised lv S a ha h) ≤ thetaLogVol (trivialised lv S a ha h) := by
  refine ⟨?_, le_of_eq (trivialised_thetaLogVol_eq_qLogVol lv S a ha h).symm⟩
  rw [trivialised_thetaLogVol_eq_qLogVol lv S a ha h, qLogVol]
  show lv S ≠ ⊤
  rw [h]
  exact WithTop.coe_ne_top

/-- 「i.e.」の形も自明に成り立つ——`C_Θ ≥ −1` は等号の場合 `C_Θ = −1` で達成される。
すなわち退化した witness は、**abc へ効く形の主張も**空虚に満たす。 -/
theorem trivialised_satisfies_CTheta (ha : 0 < a) (h : lv S = ((-a : ℝ) : WithTop ℝ))
    (C : ℝ) (hC : thetaLogVol (trivialised lv S a ha h) ≤ ((C * a : ℝ) : WithTop ℝ)) :
    -1 ≤ C := by
  rw [trivialised_thetaLogVol_eq_qLogVol lv S a ha h, qLogVol] at hC
  show (-1 : ℝ) ≤ C
  have h4 : -a ≤ C * a := by
    have : ((-a : ℝ) : WithTop ℝ) ≤ ((C * a : ℝ) : WithTop ℝ) := by
      rw [← h]; exact hC
    exact_mod_cast this
  nlinarith [ha]

-- 退化 witness であるにもかかわらず、受理ゲートを全部通ることの確認
#print axioms trivialised_satisfies_cor_3_12

/-! ## ★★ 自明化より強い事実 — 現在の `Interface` の下では `cor_3_12` は **偽** である

上の退化は「statement が無内容になりうる」ことを示した。しかし調べてみると、
事態はもっと強い。**結論の前半(有限性)を破る `PilotObjectData` が作れる**——
つまり `cor_3_12` は ∀-命題として**反証される**。

作り方: 可能な像を `{∅}` の1つだけにし(`⋃₀ {∅} = ∅`)、
対数体積を「空集合には `⊤`、それ以外には `−1`」とする。
`Interface` が課しているのは `qAbs_pos` と `qLogVol_eq` の2つだけで、
`hull` にも `possibleThetaImages` にも何の条件も無いので、これは合法な witness である。
-/

open Classical in
/-- 「空集合の対数体積は `+∞`、それ以外は `−1`」。`Interface` はこれを禁じていない。 -/
noncomputable def lvUnit (s : Set Unit) : WithTop ℝ :=
  if () ∈ s then ((-1 : ℝ) : WithTop ℝ) else ⊤

/-- ★**有限性を破る witness**。q 側は `−1`(有限)だが、Θ 側は `⊤` になる。 -/
noncomputable def finitenessCounterexample : PilotObjectData where
  Amb := Unit
  logVol := lvUnit
  hull := id
  possibleThetaImages := {∅}
  qImage := Set.univ
  qAbs := 1
  qAbs_pos := one_pos
  qLogVol_eq := by simp [lvUnit]

theorem finitenessCounterexample_thetaLogVol :
    thetaLogVol finitenessCounterexample = ⊤ := by
  simp [thetaLogVol, finitenessCounterexample, lvUnit]

/-- ★★**`cor_3_12` は現在の `Interface` の下では偽**。`sorry` は原理的に埋まらない。

これは「自明になる」より強い。PLAN §5-2 のトリアージでは**既定は①(我々のモデル化の誤り)**
——`Interface` が弱すぎて、原文が Theorem 3.11 から受け取っているはずの内容
(可能な像の和集合の正則包の対数体積が有限であること)を運べていない。

原文の有限性は Corollary 3.12 の**結論**なので、仮説として置くことはできない。
すなわち **Corollary 3.12 の前半は、Corollary 3.12 自身の逐語からは出せない**——
Theorem 3.11 の中身が要る。これが「the situation of Theorem 3.11 が未列挙」
(構造化の `<p class="open">`)の、最初の具体的な代償である。 -/
theorem cor_3_12_refutable_under_current_interface :
    ¬ ∀ D : PilotObjectData, thetaLogVol D ≠ ⊤ ∧ qLogVol D ≤ thetaLogVol D := by
  intro h
  exact (h finitenessCounterexample).1 finitenessCounterexample_thetaLogVol

#print axioms cor_3_12_refutable_under_current_interface

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
