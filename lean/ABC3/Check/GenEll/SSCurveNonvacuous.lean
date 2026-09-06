/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.SemistableFromC4
import ABC3.Found.GaloisRep.PointValuation
import ABC3.Found.GaloisRep.DegInfBaseChange
import Mathlib.RingTheory.Ideal.Int
import Mathlib.Tactic.NormNum.Prime

/-!
# 界面の測定 —— **`SSCurve` は空でない**（`Check`）

**これは原典の主張ではない**（我々の界面についての事実）ので `.src` を持たない。

## ★★★★★★★★★★★★★★★★なぜこれが要るか

`Found/GenEll/EllModuliObjects.lean` の `SSCurve`（数体 `K ⊆ ℂ` と、その上の
全ての有限素点で半安定な楕円曲線）は `EllModuliData` の `Curve` 欄の土台である。

★★☆**もし `SSCurve` が空ならば**、`Skeleton/GenEll/VeluSemistable.lean` の
`semistableAt_veluQuot_ss` も `veluQuotOK_all` も**恒真**であり、
それに依る下流（`lcyclicExc` 欄、`Lemma 3.7` の安定直線の側）は
**まとめて空虚**になる。☆2026-09-06 の測定では、木の中に
`: SSCurve :=` の形の項は `torsionExt E := E`（恒等写像）しか無かった
——すなわち witness が 1 つも無かった。

★本ファイルはその穴を塞ぐ。`Check/` 以下の `*Degenerate.lean` の 9 件が
「落とした条件は主張を偽にするか自明にする」を測っているのに対し、
本ファイルは**非空虚性の側**を測る。

## ★★★★★★★★何を作ったか —— 11a3（`y² + y = x³ − x²`）

| 量 | 値 |
|---|---|
| 定義体 | `K = ⊥ : IntermediateField ℚ ℂ`（`ℂ` の中の `ℚ`、`deg = 1`） |
| 係数 | `a₁ = 0, a₂ = −1, a₃ = 1, a₄ = 0, a₆ = 0` |
| `Δ` | `−11`（★`≠ 0` なので本当に楕円） |
| `c₄` | `16` |
| `j` | `c₄³/Δ = −4096/11` |

★★半安定性は **`gcd(Δ, c₄) = gcd(11, 16) = 1`** だけから出る:

* `v_p(11) = 0` の素点 —— `v_p(Δ) = 0` で整モデルなので `minDeltaExp p = 0`（良還元）
* `v_p(11) ≠ 0` の素点 —— Bézout `(−2)·16 + 3·11 = 1` より `v_p(c₄) = 0`、
  すなわち `c₄` が単元の整モデル（乗法還元）

★★★☆**`Δ ≠ 0` を立てるのは本質である**——`minDeltaExp` は
`if W.Δ = 0 then 0` と定義されているので、`Δ = 0` の「曲線」では
`SemistableAt` が恒真になってしまい、witness にならない。
`Δ = −11 ≠ 0` はそれを避けている。

## ★★★★★★★★★★★★おまけ —— `DegCurve` も空でない

`EllModuliData` の `Curve` 欄は `SSCurve` そのものではなく
**乗法還元を持つ**もの（`DegCurve`）に制限されている
（`Check/GenEll/EllModuliDegInfPos.lean`、第 745 の測定による）。
★11a3 は `11` の上の素点で `v_p(Δ_min) = v_p(11) > 0` なので、そこも同時に埋まる
——したがって `RealizedClass`（`EllModuliData` の対象）も空でない。

☆`11` の上の素点は going-up（`Ideal.exists_ideal_over_maximal_of_isIntegral`）で取る。

## ☆逸脱の記録

☆**無し**。本ファイルは既存の定義（`SSCurve`・`DegCurve`）をそのまま満たす項を作るだけで、
定義の側には一切触れていない。

## ☆同じ形の測定

* `Check/GenEll/NorthcottProjModelNonvacuous.lean` —— 抽象形の仮定の束が `ℙ¹` で満たされる
* `Check/GenEll/Prop17Witness.lean`・`Check/GenEll/Thm21Witness.lean`
* 本ファイル（2026-09-06）—— `SSCurve` と `DegCurve` の witness
-/

namespace ABC3.Check.GenEll

open ABC3.Found.GenEll ABC3.Found.GaloisRep
open NumberField IsDedekindDomain WeierstrassCurve IntermediateField

/-! ## ★★★★補助 —— `valAdd` の言葉での Bézout と、整モデルの良還元 -/

section Aux

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★**Bézout の関係にある 2 元は、どちらかが単元**（`valAdd` の言葉で）。

☆`Found/GenEll/PrimeOverL.lean` の `isUnit_natCast_of_coprime` は同じ内容を
`IsUnit ((n : primeSubring p))` の言葉で述べているが、そこへ橋を架けるには
`Found/GenEll/VeluBadPrimeAssembly.lean`（Vélu の組み立て一式）を import する要がある。
★本ファイルは `ValAtLeast` の超距離不等式（`Found/GaloisRep/PointValuation.lean`）だけで
同じことを言う。 -/
theorem valAdd_eq_zero_or_of_bezout (p : HeightOneSpectrum (𝓞 L)) {x y a b : L}
    (hx : x ≠ 0) (hy : y ≠ 0)
    (hxm : x ∈ primeSubring p) (hym : y ∈ primeSubring p)
    (ham : a ∈ primeSubring p) (hbm : b ∈ primeSubring p)
    (hb : a * x + b * y = 1) :
    valAdd p (Units.mk0 x hx) = 0 ∨ valAdd p (Units.mk0 y hy) = 0 := by
  rcases eq_or_ne (valAdd p (Units.mk0 x hx)) 0 with h | h1
  · exact Or.inl h
  rcases eq_or_ne (valAdd p (Units.mk0 y hy)) 0 with h | h2
  · exact Or.inr h
  exfalso
  have hx0 : 0 ≤ valAdd p (Units.mk0 x hx) := valAtLeast_of_mem hxm hx
  have hy0 : 0 ≤ valAdd p (Units.mk0 y hy) := valAtLeast_of_mem hym hy
  have hx1 : ValAtLeast p 1 x := fun _ => by omega
  have hy1 : ValAtLeast p 1 y := fun _ => by omega
  have hax : ValAtLeast p 1 (a * x) := by
    simpa using valAtLeast_mul (valAtLeast_of_mem ham) hx1
  have hby : ValAtLeast p 1 (b * y) := by
    simpa using valAtLeast_mul (valAtLeast_of_mem hbm) hy1
  have hsum : ValAtLeast p 1 (1 : L) := by
    rw [← hb]; exact valAtLeast_add hax hby
  have hfin := hsum one_ne_zero
  have hone : Units.mk0 (1 : L) one_ne_zero = 1 := Units.ext rfl
  rw [hone, valAdd_one] at hfin
  omega

/-- ★★★**整モデルで `v_p(Δ) = 0` なら `minDeltaExp p = 0`**（良還元）。

★極小モデルの判別式の付値は整モデルのそれ以下（`minDeltaExp_le_of_isIntegral`）で、
下は `0`（`minDeltaExp_nonneg`）である。 -/
theorem minDeltaExp_eq_zero_of_integral (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] (hΔ : W.Δ ≠ 0)
    (h : valAdd p (Units.mk0 W.Δ hΔ) = 0) : minDeltaExp p W = 0 := by
  have hint : WeierstrassCurve.IsIntegral (primeSubring p) ((1 : VariableChange L) • W) := by
    rwa [one_smul]
  have hle := minDeltaExp_le_of_isIntegral p W hΔ 1 hint
  have heq : valAdd p (Units.mk0 (((1 : VariableChange L) • W).Δ)
      (variableChange_Delta_ne_zero W hΔ 1)) = valAdd p (Units.mk0 W.Δ hΔ) := by
    refine valAdd_eq_of_valuation_eq p _ _ ?_
    show (p.valuation L) (((1 : VariableChange L) • W).Δ) = (p.valuation L) W.Δ
    rw [one_smul]
  rw [heq, h] at hle
  exact le_antisymm hle (minDeltaExp_nonneg p W)

/-- ★★★**極小モデルなら `minDeltaExp p = v_p(Δ)`**（変数変換が `1` の場合）。 -/
theorem minDeltaExp_eq_valAdd_of_isMinimal (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0)
    (hmin : WeierstrassCurve.IsMinimal (primeSubring p) W) :
    minDeltaExp p W = valAdd p (Units.mk0 W.Δ hΔ) := by
  have hmin' : WeierstrassCurve.IsMinimal (primeSubring p) ((1 : VariableChange L) • W) := by
    rwa [one_smul]
  rw [minDeltaExp_eq p W hΔ 1 hmin']
  refine valAdd_eq_of_valuation_eq p _ _ ?_
  show (p.valuation L) (((1 : VariableChange L) • W).Δ) = (p.valuation L) W.Δ
  rw [one_smul]

/-- ★★★**素イデアルに入る元の `valAdd` は `0` でない**。 -/
theorem valAdd_ne_zero_of_mem (p : HeightOneSpectrum (𝓞 L)) {z : 𝓞 L} (hz : z ∈ p.asIdeal)
    {x : L} (hx : x ≠ 0) (hxz : x = (z : L)) :
    valAdd p (Units.mk0 x hx) ≠ 0 := by
  intro h
  have hlt : (p.valuation L) x < 1 := by
    rw [hxz]
    exact (HeightOneSpectrum.valuation_lt_one_iff_mem p z).2 hz
  have hle : valAdd p (Units.mk0 x hx) ≤ valAdd p (1 : Lˣ) := by
    rw [valAdd_one, h]
  have h1 : (p.valuation L) ((1 : Lˣ) : L) ≤ (p.valuation L) x :=
    (valAdd_le_iff p (Units.mk0 x hx) 1).1 hle
  rw [Units.val_one, map_one] at h1
  exact absurd hlt (not_lt.2 h1)

end Aux

/-! ## ★★★★★★★★★★★★定義体 —— `ℂ` の中の `ℚ` -/

/-- ★★**定義体**——`ℂ` の中の `ℚ`（`⊥`）。 -/
noncomputable abbrev fldQ : IntermediateField ℚ ℂ := ⊥

noncomputable instance : NumberField fldQ := ⟨⟩

theorem eleven_ne_zero : (11 : fldQ) ≠ 0 := by norm_num

theorem sixteen_ne_zero : (16 : fldQ) ≠ 0 := by norm_num

/-! ## ★★★★★★★★★★★★★★★★曲線 11a3 -/

/-- ★★★★★★**`y² + y = x³ − x²`**（Cremona の `11a3`、導手 `11`）。 -/
noncomputable def curve11a3 : WeierstrassCurve fldQ :=
  { a₁ := 0, a₂ := -1, a₃ := 1, a₄ := 0, a₆ := 0 }

/-- ★★★**判別式は `−11`**。 -/
theorem curve11a3_Δ : curve11a3.Δ = -11 := by
  simp [curve11a3, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆, WeierstrassCurve.b₈]
  norm_num

/-- ★★★**`c₄ = 16`**。 -/
theorem curve11a3_c₄ : curve11a3.c₄ = 16 := by
  simp [curve11a3, WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄]
  norm_num

/-- ★★★★★**`Δ ≠ 0`**——★これが無いと `minDeltaExp` の定義（`Δ = 0` なら `0`）で
半安定性が恒真になり、witness の意味が消える。 -/
theorem curve11a3_Δ_ne_zero : curve11a3.Δ ≠ 0 := by
  rw [curve11a3_Δ]; norm_num

theorem curve11a3_c₄_ne_zero : curve11a3.c₄ ≠ 0 := by
  rw [curve11a3_c₄]; norm_num

/-- ★★★★**楕円曲線である**（`Δ` が単元）。 -/
noncomputable instance : curve11a3.IsElliptic :=
  ⟨isUnit_iff_ne_zero.2 curve11a3_Δ_ne_zero⟩

/-- ★★★**係数が整**——`0, −1, 1` はどの付値環にも入る。 -/
instance isIntegral_curve11a3 (p : HeightOneSpectrum (𝓞 fldQ)) :
    WeierstrassCurve.IsIntegral (primeSubring p) curve11a3 :=
  isIntegral_of_mem _ _ (zero_mem _) (neg_mem (one_mem _)) (one_mem _) (zero_mem _) (zero_mem _)

/-! ## ★★★★★★★★★★★★★★★★付値の 3 つの比較 -/

/-- ★★`v_p(Δ) = v_p(11)`（符号は付値に効かない）。 -/
theorem valAdd_Δ_eq (p : HeightOneSpectrum (𝓞 fldQ)) :
    valAdd p (Units.mk0 curve11a3.Δ curve11a3_Δ_ne_zero)
      = valAdd p (Units.mk0 (11 : fldQ) eleven_ne_zero) := by
  refine valAdd_eq_of_valuation_eq p _ _ ?_
  show (p.valuation fldQ) curve11a3.Δ = (p.valuation fldQ) (11 : fldQ)
  rw [curve11a3_Δ, Valuation.map_neg]

/-- ★★`v_p(c₄) = v_p(16)`。 -/
theorem valAdd_c₄_eq (p : HeightOneSpectrum (𝓞 fldQ)) :
    valAdd p (Units.mk0 curve11a3.c₄ curve11a3_c₄_ne_zero)
      = valAdd p (Units.mk0 (16 : fldQ) sixteen_ne_zero) := by
  refine valAdd_eq_of_valuation_eq p _ _ ?_
  show (p.valuation fldQ) curve11a3.c₄ = (p.valuation fldQ) (16 : fldQ)
  rw [curve11a3_c₄]

/-- ★★★★★★**`16` と `11` は同時に非単元になれない**——Bézout `(−2)·16 + 3·11 = 1`。 -/
theorem valAdd_sixteen_or_eleven (p : HeightOneSpectrum (𝓞 fldQ)) :
    valAdd p (Units.mk0 (16 : fldQ) sixteen_ne_zero) = 0 ∨
      valAdd p (Units.mk0 (11 : fldQ) eleven_ne_zero) = 0 := by
  have hm16 : (16 : fldQ) ∈ primeSubring p := by
    exact_mod_cast natCast_mem (primeSubring p) 16
  have hm11 : (11 : fldQ) ∈ primeSubring p := by
    exact_mod_cast natCast_mem (primeSubring p) 11
  have hma : (-2 : fldQ) ∈ primeSubring p := by
    exact_mod_cast intCast_mem (primeSubring p) (-2)
  have hmb : (3 : fldQ) ∈ primeSubring p := by
    exact_mod_cast natCast_mem (primeSubring p) 3
  have hbez : (-2 : fldQ) * 16 + (3 : fldQ) * 11 = 1 := by norm_num
  exact valAdd_eq_zero_or_of_bezout p sixteen_ne_zero eleven_ne_zero hm16 hm11 hma hmb hbez

/-! ## ★★★★★★★★★★★★★★★★半安定性 -/

/-- ★★★★★★★★★★★★★★★★★★★★
**11a3 はすべての有限素点で半安定**。

★機構は `gcd(Δ, c₄) = gcd(11, 16) = 1` の 1 点だけである:
`11` が単元でない素点では Bézout から `16` が単元になり、
`c₄` が単元の整モデル（＝乗法還元）として `SemistableAt` の第 2 の選言肢が立つ。
それ以外の素点では `v_p(Δ) = v_p(11) = 0` で良還元（第 1 の選言肢）である。 -/
theorem semistableAt_curve11a3 (p : HeightOneSpectrum (𝓞 fldQ)) :
    SemistableAt p curve11a3 := by
  rcases valAdd_sixteen_or_eleven p with h | h
  · -- `c₄ = 16` が単元 —— 乗法還元
    exact semistableAt_of_c4_valAdd_zero p curve11a3 curve11a3_Δ_ne_zero
      curve11a3_c₄_ne_zero ((valAdd_c₄_eq p).trans h)
  · -- `Δ = −11` が単元 —— 良還元
    exact Or.inl (minDeltaExp_eq_zero_of_integral p curve11a3 curve11a3_Δ_ne_zero
      ((valAdd_Δ_eq p).trans h))

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★`SSCurve` の witness -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`SSCurve` の具体的な元**——`ℚ` 上の 11a3。 -/
noncomputable def ssCurve11a3 : SSCurve :=
  { K := fldQ, W := curve11a3, ss := semistableAt_curve11a3 }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`SSCurve` は空でない**——★これが下流の非空虚性の錨である。 -/
theorem nonempty_SSCurve : Nonempty SSCurve := ⟨ssCurve11a3⟩

/-- ★★★**定義体の次数は `1`**（`⊥ = ℚ`）。 -/
theorem ssCurve11a3_deg : ssCurve11a3.deg = 1 := by
  show Module.finrank ℚ fldQ = 1
  exact IntermediateField.finrank_bot

/-- ★★★**`j` 不変量は `−4096/11`**（定義体の中で）。 -/
theorem curve11a3_j : curve11a3.j = -4096 / 11 := by
  have hΔ' : (curve11a3.Δ' : fldQ) = -11 := by
    rw [WeierstrassCurve.coe_Δ']; exact curve11a3_Δ
  rw [WeierstrassCurve.j, Units.val_inv_eq_inv_val, hΔ', curve11a3_c₄]
  norm_num

/-- ★★★★**`j` 不変量は `−4096/11`**（`ℂ` の中で）——★整数でないので、
この曲線は「至る所良還元」ではありえない（下の `hasMultRed` の裏付け）。 -/
theorem ssCurve11a3_j : ssCurve11a3.j = -4096 / 11 := by
  have h1 : ((4096 : fldQ) : ℂ) = 4096 := map_ofNat (fldQ.val.toRingHom) 4096
  have h2 : ((11 : fldQ) : ℂ) = 11 := map_ofNat (fldQ.val.toRingHom) 11
  show ((curve11a3.j : fldQ) : ℂ) = -4096 / 11
  rw [curve11a3_j]
  push_cast
  rw [h1, h2]

/-! ## ★★★★★★★★★★★★★★`11` の上の素点 -/

/-- ★★★★★**`11` の上に素点がある**——going-up。 -/
theorem exists_prime_over_eleven :
    ∃ p : HeightOneSpectrum (𝓞 fldQ), ((11 : ℤ) : 𝓞 fldQ) ∈ p.asIdeal := by
  have hinj : Function.Injective (algebraMap ℤ (𝓞 fldQ)) :=
    RingHom.injective_int (algebraMap ℤ (𝓞 fldQ))
  have hker : RingHom.ker (algebraMap ℤ (𝓞 fldQ)) ≤ Ideal.span {(11 : ℤ)} := by
    rw [(RingHom.injective_iff_ker_eq_bot _).1 hinj]
    exact bot_le
  haveI : Fact (Nat.Prime 11) := ⟨by norm_num⟩
  haveI hmax : (Ideal.span {(11 : ℤ)}).IsMaximal := by
    have h := Int.ideal_span_isMaximal_of_prime 11
    simpa using h
  obtain ⟨Q, hQmax, hQcomap⟩ :=
    Ideal.exists_ideal_over_maximal_of_isIntegral (Ideal.span {(11 : ℤ)}) hker
  have hmem : ((11 : ℤ) : 𝓞 fldQ) ∈ Q := by
    have h11 : (11 : ℤ) ∈ Ideal.span {(11 : ℤ)} := Ideal.mem_span_singleton_self _
    rw [← hQcomap] at h11
    exact h11
  have hne : Q ≠ ⊥ := by
    intro hQ
    rw [hQ, ← RingHom.ker_eq_comap_bot,
      (RingHom.injective_iff_ker_eq_bot _).1 hinj] at hQcomap
    have h0 : (11 : ℤ) ∈ (⊥ : Ideal ℤ) := by
      rw [hQcomap]; exact Ideal.mem_span_singleton_self _
    rw [Ideal.mem_bot] at h0
    norm_num at h0
  exact ⟨⟨Q, hQmax.isPrime, hne⟩, hmem⟩

/-! ## ★★★★★★★★★★★★★★★★★★`DegCurve` の witness -/

/-- ★★★★★★★★★★★★★★★★★★★★
**11a3 は乗法還元を持つ**——`11` の上の素点で `v_p(Δ_min) = v_p(11) > 0`。

★その素点では `c₄ = 16` が単元（Bézout）なので与えたモデルが極小であり、
`minDeltaExp p = v_p(Δ) ≠ 0` である。 -/
theorem hasMultRed_ssCurve11a3 : ssCurve11a3.HasMultRed := by
  obtain ⟨p, hp⟩ := exists_prime_over_eleven
  refine ⟨p, ?_⟩
  show minDeltaExp p curve11a3 ≠ 0
  -- `11` は `p` で単元でない
  have h11 : valAdd p (Units.mk0 (11 : fldQ) eleven_ne_zero) ≠ 0 := by
    refine valAdd_ne_zero_of_mem p hp eleven_ne_zero ?_
    have h1 : ((11 : ℤ) : 𝓞 fldQ) = (11 : 𝓞 fldQ) := by norm_num
    rw [h1]
    exact (map_ofNat (algebraMap (𝓞 fldQ) fldQ) 11).symm
  -- したがって `c₄ = 16` は単元、すなわちモデルは極小
  have h16 : valAdd p (Units.mk0 (16 : fldQ) sixteen_ne_zero) = 0 := by
    rcases valAdd_sixteen_or_eleven p with h | h
    · exact h
    · exact absurd h h11
  have hmin : WeierstrassCurve.IsMinimal (primeSubring p) curve11a3 :=
    isMinimal_of_c4_valAdd_eq_zero p curve11a3 curve11a3_Δ_ne_zero curve11a3_c₄_ne_zero
      ((valAdd_c₄_eq p).trans h16)
  rw [minDeltaExp_eq_valAdd_of_isMinimal p curve11a3 curve11a3_Δ_ne_zero hmin, valAdd_Δ_eq p]
  exact h11

/-- ★★★★★★★★★★★★★★★★**`DegCurve` の具体的な元**。 -/
noncomputable def degCurve11a3 : DegCurve :=
  { toSSCurve := ssCurve11a3, multRed := hasMultRed_ssCurve11a3 }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`DegCurve` は空でない**——★`EllModuliData` の `Curve` 欄（乗法還元を持つ半安定曲線）が
実際に対象を持つ。 -/
theorem nonempty_DegCurve : Nonempty DegCurve := ⟨degCurve11a3⟩

/-- ★★★★★★★★★★★★★★★★★★★★
**`RealizedClass` も空でない**——`EllModuliData` の対象の全体が空でない。 -/
theorem nonempty_RealizedClass : Nonempty RealizedClass :=
  ⟨⟨degCurve11a3.j, ⟨degCurve11a3, rfl⟩⟩⟩

/-- ★★★★★★★★★★★★**`deg∞ > 0` の類が実際に在る**——
`Check/GenEll/EllModuliDegInfPos.lean`（第 745）が界面から強いた条件を、
具体的な曲線が満たしている。 -/
theorem degInfJ_pos_ssCurve11a3 : 0 < degInfJ ssCurve11a3.j :=
  degCurve11a3.degInf_pos

/-!
## ★★★★★★★★★★★★★★★★★★★★★★★★★★
界面の測定（続き、2026-09-06）—— `VeluQuotOK` の内側の ∀ も空虚でない

★上で `SSCurve` が空でないことは確定したが、`Skeleton/GenEll/VeluSemistable.lean` の
`veluQuotOK_all` が主張する `VeluQuotOK E l`
（`Found/GenEll/Lemma37Hdag/Lemma35.lean`）は

  ∀ (M : IntermediateField E.fld E.alg) [FiniteDimensional E.fld M]
    (Q : (E.W.baseChange M).toAffine.Point), addOrderOf Q = l → …

という形の ∀ である。★★☆したがって「そのような `M` と `Q` が 1 つも無い」なら
`VeluQuotOK` は前提が空で恒真になり、`Lemma 3.5` の内容が消える。
本節はその側を塞ぐ——11a3 の有理 5-捩れ点 `(0, 0)` を実際に構成する。

### ★★★★★★★★どう計算したか（`l = 5` の witness）

`y² + y = x³ − x²` の上で `negY x y = −y − 1` であり、`P = (0, 0)` から

| 倍 | 座標 | 使った公式 |
|---|---|---|
| `P` | `(0, 0)` | `nonsingular_zero`（`a₆ = 0` と `a₃ = 1 ≠ 0`） |
| `2P` | `(1, −1)` | 接線の傾き `λ = (3x² + 2a₂x + a₄ − a₁y)/(2y + a₁x + a₃) = 0/1 = 0` |
| `4P` | `(0, −1)` | 接線の傾き `λ = (3 − 2)/(−2 + 1) = −1` |
| `5P` | `O` | `4P` と `P` は `x` が同じで `y₁ = negY x₂ y₂`（縦線） |

★`addOrderOf = 5` は `addOrderOf_eq_prime`（`5` は素数・`5 • P = 0`・`P ≠ 0`）で出る。
すなわち `2P, 3P, 4P ≠ 0` を別に計算する要は無い。

### ★★★★★★★★なぜ体を一般化したか

`VeluQuotOK` の `Q` は `E.W.baseChange M` の上の点なので、
`fldQ` の上だけで計算しても届かない。★上の表の計算は「`a₁ = 0, a₂ = −1, a₃ = 1,
a₄ = 0, a₆ = 0` という係数だけ」を使い、体の性質は `1 ≠ 0` しか使わない。
そこで `model11a3 L`（任意の体 `L` の上の同じ Weierstrass 方程式）で計算し、
`L = fldQ` と `L = M` の両方に流す。
☆`(model11a3 L).baseChange A = model11a3 A` は係数が `0, ±1` だけなので
`map_zero`・`map_one`・`map_neg` で出る。

### ☆逸脱の記録

☆無し。既存の `curve11a3` は書き換えていない
（`curve11a3 = model11a3 fldQ` は `rfl` で、両者は同じ項である）。
-/

section Model

variable (L : Type) [Field L]

/-- ★★★★★★**`y² + y = x³ − x²`**——任意の体の上の同じ Weierstrass 方程式。

☆`curve11a3` はこれの `L = fldQ` の場合そのものである（`curve11a3_eq_model11a3`）。 -/
def model11a3 : WeierstrassCurve L :=
  { a₁ := 0, a₂ := -1, a₃ := 1, a₄ := 0, a₆ := 0 }

theorem model11a3_a₁ : (model11a3 L).toAffine.a₁ = 0 := rfl

theorem model11a3_a₂ : (model11a3 L).toAffine.a₂ = -1 := rfl

theorem model11a3_a₃ : (model11a3 L).toAffine.a₃ = 1 := rfl

theorem model11a3_a₄ : (model11a3 L).toAffine.a₄ = 0 := rfl

theorem model11a3_a₆ : (model11a3 L).toAffine.a₆ = 0 := rfl

/-- ★★★**`−(x, y) = (x, −y − 1)`**（`a₁ = 0`・`a₃ = 1` なので）。 -/
theorem model11a3_negY (x y : L) : (model11a3 L).toAffine.negY x y = -y - 1 := by
  rw [WeierstrassCurve.Affine.negY, model11a3_a₁, model11a3_a₃]; ring

/-- ★★`(0, 0)` は 2-捩れでない（接線が縦でない）。 -/
theorem model11a3_Y_ne_negY_zero : (0 : L) ≠ (model11a3 L).toAffine.negY 0 0 := by
  rw [model11a3_negY]
  intro h
  exact one_ne_zero (α := L) (by linear_combination h)

/-- ★★`(1, −1)` は 2-捩れでない。 -/
theorem model11a3_Y_ne_negY_one : (-1 : L) ≠ (model11a3 L).toAffine.negY 1 (-1) := by
  rw [model11a3_negY]
  intro h
  exact one_ne_zero (α := L) (by linear_combination -h)

/-- ★★★★**`(0, 0)` は曲線上の非特異点**——`a₆ = 0` と `a₃ = 1 ≠ 0`。 -/
theorem nonsingular_model11a3_zero_zero : (model11a3 L).toAffine.Nonsingular 0 0 := by
  rw [WeierstrassCurve.Affine.nonsingular_zero, model11a3_a₃, model11a3_a₆]
  exact ⟨rfl, Or.inl one_ne_zero⟩

/-- ★★★**`(1, −1)` は曲線上の非特異点**（`= 2P`）。 -/
theorem nonsingular_model11a3_one_negOne : (model11a3 L).toAffine.Nonsingular 1 (-1) := by
  rw [WeierstrassCurve.Affine.nonsingular_iff, WeierstrassCurve.Affine.equation_iff,
    model11a3_a₁, model11a3_a₂, model11a3_a₃, model11a3_a₄, model11a3_a₆]
  exact ⟨by ring, Or.inr fun h => one_ne_zero (α := L) (by linear_combination -h)⟩

/-- ★★★**`(0, −1)` は曲線上の非特異点**（`= −P = 4P`）。 -/
theorem nonsingular_model11a3_zero_negOne : (model11a3 L).toAffine.Nonsingular 0 (-1) := by
  have h := (WeierstrassCurve.Affine.nonsingular_neg
    (W' := (model11a3 L).toAffine) 0 0).2 (nonsingular_model11a3_zero_zero L)
  rw [model11a3_negY] at h
  rwa [show (-0 - 1 : L) = -1 by ring] at h

/-- ★★**座標が一致すれば点も一致**（非特異性の証明は `Prop` なので効かない）。 -/
theorem model11a3_point_ext {x y x' y' : L} (h : (model11a3 L).toAffine.Nonsingular x y)
    (h' : (model11a3 L).toAffine.Nonsingular x' y') (hx : x = x') (hy : y = y') :
    WeierstrassCurve.Affine.Point.some x y h = WeierstrassCurve.Affine.Point.some x' y' h' := by
  subst hx; subst hy; rfl

/-- ★★★★★★**`P = (0, 0)`**——11a3 の有理 5-捩れ点。 -/
def model11a3P : (model11a3 L).toAffine.Point :=
  WeierstrassCurve.Affine.Point.some 0 0 (nonsingular_model11a3_zero_zero L)

/-- ★★★**`2P = (1, −1)`**。 -/
def model11a3P2 : (model11a3 L).toAffine.Point :=
  WeierstrassCurve.Affine.Point.some 1 (-1) (nonsingular_model11a3_one_negOne L)

/-- ★★★**`4P = (0, −1)`**（`= −P`）。 -/
def model11a3P4 : (model11a3 L).toAffine.Point :=
  WeierstrassCurve.Affine.Point.some 0 (-1) (nonsingular_model11a3_zero_negOne L)

/-- ★★★★**`P ≠ O`**——★これが `addOrderOf` を `1` でなく `5` にする。 -/
theorem model11a3P_ne_zero : model11a3P L ≠ 0 :=
  WeierstrassCurve.Affine.Point.some_ne_zero (nonsingular_model11a3_zero_zero L)

/-- ★★★★★**係数が `0, ±1` だけなので、底の変換で式は変わらない**。 -/
theorem model11a3_baseChange (A : Type) [Field A] [Algebra L A] :
    (model11a3 L).baseChange A = model11a3 A := by
  ext <;> simp [model11a3, WeierstrassCurve.baseChange, WeierstrassCurve.map]

variable [DecidableEq L]

/-- ★★★**`P` での接線の傾きは `0`**——`(3·0² + 2(−1)·0 + 0 − 0·0)/(0 − (−0 − 1)) = 0/1`。 -/
theorem model11a3_slope_zero : (model11a3 L).toAffine.slope 0 0 0 0 = 0 := by
  rw [WeierstrassCurve.Affine.slope_of_Y_ne rfl (model11a3_Y_ne_negY_zero L),
    model11a3_negY, model11a3_a₁, model11a3_a₂, model11a3_a₄,
    show ((3 : L) * 0 ^ 2 + 2 * (-1) * 0 + 0 - 0 * 0) = 0 by ring,
    show ((0 : L) - (-0 - 1)) = 1 by ring, div_one]

/-- ★★★**`2P` での接線の傾きは `−1`**——`(3 − 2)/(−1) = −1`。 -/
theorem model11a3_slope_one : (model11a3 L).toAffine.slope 1 1 (-1) (-1) = -1 := by
  rw [WeierstrassCurve.Affine.slope_of_Y_ne rfl (model11a3_Y_ne_negY_one L),
    model11a3_negY, model11a3_a₁, model11a3_a₂, model11a3_a₄,
    show ((3 : L) * 1 ^ 2 + 2 * (-1) * 1 + 0 - 0 * (-1)) = 1 by ring,
    show ((-1 : L) - (-(-1) - 1)) = -1 by ring, div_neg, div_one]

/-- ★★★★★**`P + P = 2P`**——`addX 0 0 0 = 1`, `addY 0 0 0 0 = −1`。 -/
theorem model11a3_add_P_P : model11a3P L + model11a3P L = model11a3P2 L := by
  rw [model11a3P, model11a3P2,
    WeierstrassCurve.Affine.Point.add_self_of_Y_ne (model11a3_Y_ne_negY_zero L)]
  refine model11a3_point_ext L _ _ ?_ ?_
  · rw [WeierstrassCurve.Affine.addX, model11a3_slope_zero, model11a3_a₁, model11a3_a₂]; ring
  · rw [WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negAddY,
      WeierstrassCurve.Affine.addX, model11a3_slope_zero, model11a3_negY,
      model11a3_a₁, model11a3_a₂]
    ring

/-- ★★★★★**`2P + 2P = 4P`**——`addX 1 1 (−1) = 0`, `addY 1 1 (−1) (−1) = −1`。 -/
theorem model11a3_add_P2_P2 : model11a3P2 L + model11a3P2 L = model11a3P4 L := by
  rw [model11a3P2, model11a3P4,
    WeierstrassCurve.Affine.Point.add_self_of_Y_ne (model11a3_Y_ne_negY_one L)]
  refine model11a3_point_ext L _ _ ?_ ?_
  · rw [WeierstrassCurve.Affine.addX, model11a3_slope_one, model11a3_a₁, model11a3_a₂]; ring
  · rw [WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negAddY,
      WeierstrassCurve.Affine.addX, model11a3_slope_one, model11a3_negY,
      model11a3_a₁, model11a3_a₂]
    ring

/-- ★★★★★**`4P + P = O`**——`x` が同じで `y₁ = negY x₂ y₂`（縦線）。 -/
theorem model11a3_add_P4_P : model11a3P4 L + model11a3P L = 0 := by
  rw [model11a3P4, model11a3P]
  refine WeierstrassCurve.Affine.Point.add_of_Y_eq rfl ?_
  rw [model11a3_negY]; ring

/-- ★★★★★★**`5P = O`**。 -/
theorem model11a3_five_nsmul : (5 : ℕ) • model11a3P L = 0 := by
  have h5 : (5 : ℕ) • model11a3P L
      = model11a3P L + model11a3P L + (model11a3P L + model11a3P L) + model11a3P L := by abel
  rw [h5, model11a3_add_P_P, model11a3_add_P2_P2, model11a3_add_P4_P]

/-- ★★★★★★★★★★★★★★★★
**`(0, 0)` の位数はちょうど `5`**——任意の体の上で。

★`5` は素数なので `5P = O` と `P ≠ O` だけで位数が定まる
（`addOrderOf_eq_prime`）。`2P, 3P, 4P ≠ O` を別に計算する要は無い。 -/
theorem model11a3_addOrderOf : addOrderOf (model11a3P L) = 5 := by
  haveI : Fact (Nat.Prime 5) := ⟨by norm_num⟩
  exact addOrderOf_eq_prime (model11a3_five_nsmul L) (model11a3P_ne_zero L)

/-- ★★★★★★**この方程式に等しい曲線には位数 `5` の点がある**。

☆`W` を変数に置いてから `subst` するのは、`(model11a3 L).toAffine.Point` と
`W.toAffine.Point` の間の輸送を型検査器にやらせるためである
（`W` が変数でないと `subst` が効かない）。 -/
theorem exists_point_addOrderOf_five (W : WeierstrassCurve L) (hW : W = model11a3 L) :
    ∃ Q : W.toAffine.Point, addOrderOf Q = 5 := by
  subst hW
  exact ⟨model11a3P L, model11a3_addOrderOf L⟩

end Model

/-! ## ★★★★★★★★★★★★★★★★★★★★11a3 の 5-捩れ点 -/

/-- ★★★**`curve11a3` は `model11a3 fldQ` そのもの**（同じ項）。 -/
theorem curve11a3_eq_model11a3 : curve11a3 = model11a3 fldQ := rfl

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**11a3 の有理点 `(0, 0)`**——★これが `VeluQuotOK` の `l = 5` の witness の芽である。 -/
noncomputable def point11a3 : curve11a3.toAffine.Point := model11a3P fldQ

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`(0, 0)` の位数は `5`**——古典的事実「11a3 の有理捩れ群は `ℤ/5ℤ`」の、
`5` が位数であるという側を機械にかけたもの。 -/
theorem addOrderOf_point11a3 : addOrderOf point11a3 = 5 := model11a3_addOrderOf fldQ

/-- ★★★**`(0, 0) ≠ O`**。 -/
theorem point11a3_ne_zero : point11a3 ≠ 0 := model11a3P_ne_zero fldQ

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★
`VeluQuotOK` の内側の ∀ は空虚でない -/

/-- ★★★★★★★★★★★★★★★★★★★★★★
**`ssCurve11a3` をどの拡大体へ底変換しても位数 `5` の点がある**。

★`VeluQuotOK E l` の内側の `∀ Q, addOrderOf Q = l → …` は
`E.W.baseChange M` の上の点についての ∀ である。`l = 5` のとき、
本補題はその前提を満たす `Q` が「どの `M` でも」存在することを言う。 -/
theorem exists_addOrderOf_eq_five_baseChange (M : Type) [Field M]
    [Algebra ssCurve11a3.fld M] [DecidableEq M] :
    ∃ Q : (ssCurve11a3.W.baseChange M).toAffine.Point, addOrderOf Q = 5 :=
  exists_point_addOrderOf_five M _ (model11a3_baseChange ssCurve11a3.fld M)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`VeluQuotOK ssCurve11a3 5` の内側の ∀ は前提が空でない**。

☆`M` と `DecidableEq` の取り方は `Found/GenEll/Lemma37Hdag/Lemma35.lean` の
`VeluQuotOK` の定義に合わせてある（そこでは `DecidableEq` を
`Classical.propDecidable` で入れている）。
☆`ssCurve11a3.alg` は定義上 `AlgebraicClosure ssCurve11a3.fld` なので、
本ファイルは `EllModuliGalois/Theorem38.lean` を import せずに同じ型を書いている。

★★☆いちばん安い witness は自明な拡大 `M = ⊥` である。 -/
theorem exists_veluQuotOK_premise_witness :
    ∃ (M : IntermediateField ssCurve11a3.fld (AlgebraicClosure ssCurve11a3.fld))
      (_ : FiniteDimensional ssCurve11a3.fld M),
      letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
      ∃ Q : (ssCurve11a3.W.baseChange M).toAffine.Point, addOrderOf Q = 5 := by
  refine ⟨⊥, inferInstance, ?_⟩
  exact @exists_point_addOrderOf_five _ _ (fun a b => Classical.propDecidable (a = b)) _
    (model11a3_baseChange fldQ _)

end ABC3.Check.GenEll
