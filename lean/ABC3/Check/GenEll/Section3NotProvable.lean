import ABC3.Interface.GenEll.TateLocal
import ABC3.Interface.GenEll.EllModuli
import ABC3.Found.GenEll.BDClass

/-!
# 退化封じの検査 —— **§3 の 5 つの `sorry` は「難しい」のではなく偽である**(`Check`)

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。

## ★★★★★★★★2026-08-26 の診断

`Skeleton/GenEll/Section3.lean` に残る 5 つの `sorry` は

| 定理 | 量化している界面 |
|---|---|
| `lemma_3_2` | `∀ D : TateLocalData` |
| `potLocalHeight_indep` | `∀ D : TateLocalData` ★**2026-08-26 に閉じた** |
| `prop_3_4` | `∀ D : EllModuliData` ★**2026-08-26 に閉じた**(第 361) |
| `lemma_3_5` | `∀ D : EllModuliData` |
| `lemma_3_7` | `∀ D : EllModuliData` |

という形をしている。★★**ところがこの 2 つの界面は公理を 1 つも持たない**
(`TateLocalData` は `vq_pos`・`ramIdx_pos` の 2 本、`EllModuliData` は 0 本)。
★★★したがって**退化した `D` を作れば主張は破れる**——`sorry` は
「まだ証明していない」ではなく「**この形では証明できない**」ことの印である。

## ★★★★★★これは G6 が 2026-08-17 に通った道と同じである

`Interface/GaloisRep/Reduction.lean` の `TateCurveData` は、以前
`LocalField : Type` / `Curve : LocalField → Type` と**世界ごと posit** していて
`PUnit` と定数 `1` で埋まった。★そこで mathlib の `WeierstrassCurve` と
正規化付値に**接地**し、**Tate 一意化そのもの**を要求する形に作り直した。

★★`Interface/GenEll/TateLocal.lean` は**その捨てたはずの形のまま**である
(`LocalField : Type` / `Curve : LocalField → Type`)。
★★★`EllModuliData` も `EllClass : Type` / `faltingsHeight : EllClass → ℝ` を
条件なしで持つだけである。

## ★★★★本ファイルが示すこと(限界の明示)

★示すのは「**現在の界面の上では §3 の 5 つは偽である**」ことだけである。
★★「原文の `Lemma 3.2` / `Proposition 3.4` が偽」ということでは**まったくない**。
★★★塞ぎ方も分かっている——G6 と同じく、**既に達成済みの
`TateCurveData`(G6)・`FaltingsHeightData`(G8)に接地する**ことである。

## ★★界面の欠陥の一覧(本ファイルで 2 つ増える)

| # | 場所 | 欠陥の型 | 塞いだ |
|---|---|---|---|
| 1-5 | G6/G7/G8 | 充足不能・弱すぎる | 第 302-318 |
| 6 | G8 `htFalt` | **弱すぎる**(`deg∞/12` で満たせる) | 第 357 |
| 7 | **`TateLocalData`** | **弱すぎる**(世界ごと posit、`Unit` で埋まる) | ★**部分**(第 360: `Remark 3.3.1`) |
| 8 | **`EllModuliData`** | **弱すぎる**(高さの間に関係が無い) | ★**部分**(第 361: `Proposition 3.4`) |
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll ABC3.Found.GenEll

/-! ## ★★★★★★★★1. `TateLocalData` は 2026-08-26 に接地された

★以前はここに `degenerateTateLocal`･`lemma_3_2_i_false`･`lemma_3_2_ii_false` があった
——世界を `Unit`、`v_K(q_E)` を定数 `1`、`deg_∞` を定数にして破れていた。

★★次の 5 つを欄に出したので、**その退化 witness はもはや作れない**:
`vq_baseChange`(第 360)･`logResidueCard`･`logResidueCard_pos`･
`degInf_eq`･`vq_quotMu`･`stableLine_dvd_or_cyclotomic`(第 363)。

★★★**充足不能になっていないかを確かめる**のが下の `tateLocalSatisfiable` である
——欄を強めるときは**強めすぎて充足不能にする**のが欠陥 #1･#2･#4･#5 であった。 -/

/-- ★★★★★**強めた `TateLocalData` は充足可能である**。

`v_K(q_E) = n+1`･`deg_∞ = n+1`･`E/μ_l` を `l(n+1)-1` とすればすべての欄が成り立つ。
★`l` を素数に限ったのがここで効く——`l = 0` なら `vq_quotMu` と `vq_pos` が衝突する。 -/
def tateLocalSatisfiable : TateLocalData where
  LocalField := Unit
  Curve := fun _ => ℕ
  vq := fun _ n => n + 1
  vq_pos := fun _ n => Nat.succ_pos n
  degInf := fun _ n => ((n + 1 : ℕ) : ℝ)
  StableLine := fun _ _ _ => Empty
  IsCyclotomic := fun _ => True
  quotMu := fun _ n l => l * (n + 1) - 1
  MultExt := fun _ _ => Unit
  extField := fun _ => ()
  baseChange := fun _ => 0
  ramIdx := fun _ => 1
  ramIdx_pos := fun _ => Nat.one_pos
  vq_baseChange := fun _ _ => rfl
  logResidueCard := fun _ => 1
  logResidueCard_pos := fun _ => one_pos
  degInf_eq := by
    intro K E
    show ((E + 1 : ℕ) : ℝ) = ((E + 1 : ℕ) : ℝ) * 1
    ring
  vq_quotMu := by
    intro K E l hl
    show (l * (E + 1) - 1) + 1 = l * (E + 1)
    have h1 : 1 ≤ l * (E + 1) :=
      Nat.one_le_iff_ne_zero.2 (Nat.mul_ne_zero hl.ne_zero (Nat.succ_ne_zero E))
    omega
  stableLine_dvd_or_cyclotomic := by
    intro K E l hl N
    exact N.elim

/-! ## ★★★★★★★★2. `Proposition 3.4` も 2026-08-26 に閉じた

★以前はここに `prop_3_4_first_false`･`prop_3_4_finite_false` があった——
`deg_∞ ≡ 0`･`ht_∞ = n`･`ht^Falt ≡ 0` とすれば破れていた。

★★`Interface/GenEll/EllModuli.lean` に、原文が `Proposition 3.4` の証明で
**実際に引いている 4 つ**(`degInf_le_htInf`･`htInf_bdeq_faltings`･
`faltingsHeight_bddBelow`･`northcott`)を欄として足したので、
**その退化 witness はもはや作れない**。 -/

/-! ## ★★★★★★★3. `≲` を印字どおりに読むと `ht_∞` が上に有界になる -/

/-- ★★★★★★★**逐語どおりの `BDle` で `Proposition 3.4` の鎖を書くと、
`ht_∞` が上に有界になってしまう**。

`BDle α β` は `∃C, ∀x, β x − α x ≤ C`(原典 `Definition 1.2, (ii)` の印字どおり)。
★2 番目 `BDle ht_∞ (12(1+ε)ht^Falt)` と 3 番目 `BDle (12(1+ε)ht^Falt) ((1+ε)ht_∞)` を
足すと `(1+ε)ht_∞ ≤ ht_∞ + C`、すなわち `ε·ht_∞ ≤ C` が出る。

★★モジュライ上の高さは上に有界でないので、**その読みは意図されたものではない**。
★★★`Skeleton/GenEll/Section3.lean` の `prop_3_4` は `BDge` で書き直した(逸脱を記録)。 -/
theorem bdle_chain_forces_bounded {F : Type*} (htInf faltings : F → ℝ) (eps : ℝ) (heps : 0 < eps)
    (h2 : BDle htInf (fun x => 12 * (1 + eps) * faltings x))
    (h3 : BDle (fun x => 12 * (1 + eps) * faltings x) (fun x => (1 + eps) * htInf x)) :
    ∃ C : ℝ, ∀ x, htInf x ≤ C := by
  obtain ⟨C₂, hC₂⟩ := h2
  obtain ⟨C₃, hC₃⟩ := h3
  refine ⟨(C₂ + C₃) / eps, fun x => ?_⟩
  have h2x := hC₂ x
  have h3x := hC₃ x
  rw [le_div_iff₀ heps]
  nlinarith

/-! ## ★★★★★4. `ε` 入りの 2 本は posit ではなく導出である -/

/-- ★★★★★**`ht_∞ ≈ 12·ht^Falt` と `ht^Falt` の下界から `ε` 入りの不等式が出る**。

★posit しているのは `≈` の方であって、`ε` 入りの形ではない。 -/
theorem eps_ineq_is_derived {F : Type*} (htInf faltings : F → ℝ) (eps : ℝ) (heps : 0 < eps)
    (M B : ℝ) (hM : ∀ x, |htInf x - 12 * faltings x| ≤ M) (hB : ∀ x, B ≤ faltings x) :
    BDge htInf (fun x => 12 * (1 + eps) * faltings x) := by
  refine ⟨M - 12 * eps * B, fun x => ?_⟩
  have h1 := abs_le.1 (hM x)
  have h2 := hB x
  nlinarith [h1.1, h1.2]

end ABC3.Check.GenEll
