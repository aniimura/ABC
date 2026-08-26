import ABC3.Found.GaloisRep.DetWiring

/-!
# Galois (G5) 第 207 ブロック —— **★★★★★★★★★`det ρ = 円分指標`(生成性を仮定した形)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★配線が閉じた

第 204(円分指標の段)・第 205(射影と `ℤ_l` 作用)・第 206(行列成分の読み替え)を繋ぎ、
**`det_galRep_eq_cyclotomic` を `hgen` 1 つだけを仮定して証明した**:

    σ ζ = ζ^{(toZModPow n (det (galRep W l e σ))).val}

★`hgen` は「`ζ^{l^n} = 1` なら `ζ = e_{l^n}(P,Q)^k`」——すなわち
**`e_{l^n}(P,Q)` が原始 `l^n` 乗根であること**(非退化性 + 基底性)である。

### ★★★★★機構

| 段 | 使うもの |
|---|---|
| `P := tateVec(single₀)`・`Q := tateVec(single₁)` | `T_l E` の基底の第 `n` 成分 |
| `σP = a·P + c·Q`・`σQ = b·P + d·Q` | 第 205(射影と Galois の可換)+ 第 206(列を基底で読む) |
| `σ ζ · ζ^{bc} = ζ^{ad}` | 第 204 |
| `det + bc ≡ ad (mod l^n)` | 第 206 |
| `σ ζ = ζ^{det}` | 本ブロックの `eq_pow_of_mul_eq` |

### ★★★これで (G5) は単一の入力に帰着した

`hgen` = 非退化性 + 基底性。★非退化性は第 197 で `F(E)^{E[n]} = [n]^*F(E)` に絞れ、
第 196・198 と合わせて `x(nP) = Φ_n/ΨSq_n`、すなわち
**`normEDS` が楕円列であること**(mathlib 自身の TODO)に帰着している。
★★基底性は `T_l E ≅ ℤ_l²` の基底の第 `n` 成分が `E[l^n]` の基底になること
(射影の全射性は第 186 で済み)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `eq_pow_of_mul_eq` | ★★★★★最後の打ち消し |
| `tateVec` ほか | `T_l E` の基底から `E[l^n]` への加法写像 |
| `galPoint_tateVec` | ★★★★`σ` は行列で動く |
| `det_cyclotomic_of_gen` | ★★★★★★★★★**`det ρ = 円分指標`(`hgen` 仮定)** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- ★★★★★**最後の打ち消し**——`w·ζ^{bc} = ζ^{ad}` と `D + bc ≡ ad` から `w = ζ^D`。 -/
theorem eq_pow_of_mul_eq {F : Type} [Field F] {m : ℕ} (hm : 1 ≤ m) {ζ : F} (hζ : ζ ^ m = 1)
    {w : F} {D bc ad : ℕ} (hcong : (D + bc) % m = ad % m)
    (h : w * ζ ^ bc = ζ ^ ad) : w = ζ ^ D := by
  have hζ0 : ζ ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega : m ≠ 0)] at hζ
    exact zero_ne_one hζ
  refine mul_right_cancel₀ (pow_ne_zero bc hζ0) ?_
  rw [h, ← pow_add]
  exact (pow_eq_pow_of_modEq hζ hcong).symm

/-- `T_l E` のベクトルから `E[l^n]` への加法写像。 -/
noncomputable def tateVec (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (n : ℕ) :
    (Fin 2 → ℤ_[l]) →+ (W.baseChange L).toAffine.Point :=
  ((AddSubgroup.subtype (torsionPoints (W.baseChange L) (l ^ n))).comp
    (tateProj (W.baseChange L) l n)).comp e.symm.toAddMonoidHom

theorem tateVec_apply (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (n : ℕ) (v : Fin 2 → ℤ_[l]) :
    tateVec W l e n v
      = ((tateProj (W.baseChange L) l n (e.symm v) : (W.baseChange L).toAffine.Point)) := rfl

theorem tateVec_torsion (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (n : ℕ) (v : Fin 2 → ℤ_[l]) :
    (l ^ n) • tateVec W l e n v = 0 :=
  (tateProj (W.baseChange L) l n (e.symm v)).2

/-- ★★★★**`σ` は行列で動く**。 -/
theorem galPoint_tateVec (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (n : ℕ) (σ : L ≃ₐ[K] L)
    (v : Fin 2 → ℤ_[l]) :
    galPoint W σ (tateVec W l e n v)
      = tateVec W l e n (Matrix.mulVec (galMat W l e σ) v) := by
  have hg : galTate W l σ (e.symm v) = e.symm (Matrix.mulVec (galMat W l e σ) v) := by
    refine e.injective ?_
    rw [AddEquiv.apply_symm_apply, galMat_apply, AddEquiv.apply_symm_apply]
  rw [tateVec_apply, tateVec_apply, ← hg]
  exact (tateProj_galTate W l n σ (e.symm v)).symm

/-- ★★★★★★★★★**`det ρ = 円分指標`(生成性を仮定した形)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`hgen` は「`e_{l^n}(P,Q)` が原始 `l^n` 乗根であること」——非退化性 + 基底性である。
★★それ以外はすべて第 204-206 で済んでいる。 -/
theorem det_cyclotomic_of_gen [IsAlgClosed L] [CharZero L]
    (W : WeierstrassCurve K) [((W.baseChange L).toAffine).IsElliptic]
    (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L) (n : ℕ) (ζ : L) (hζ : ζ ^ (l ^ n) = 1)
    (hgen : ∀ z : L, z ^ (l ^ n) = 1 →
      ∃ k : ℕ, (weilPairingVal (W.baseChange L).toAffine (l ^ n)
        (tateVec W l e n (Pi.single 0 1)) (tateVec W l e n (Pi.single 1 1))) ^ k = z) :
    σ ζ = ζ ^ ((PadicInt.toZModPow n
      ((galRep W l e σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).det)).val) := by
  have hm : 1 ≤ l ^ n := Nat.one_le_pow _ _ (Nat.pos_of_ne_zero (Fact.out : l.Prime).ne_zero)
  set M := galMat W l e σ with hM
  set P := tateVec W l e n (Pi.single 0 1) with hPdef
  set Q := tateVec W l e n (Pi.single 1 1) with hQdef
  have hPt : (l ^ n) • P = 0 := tateVec_torsion W l e n _
  have hQt : (l ^ n) • Q = 0 := tateVec_torsion W l e n _
  have hσP : galPoint W σ P
      = ((PadicInt.toZModPow n (M 0 0)).val) • P
        + ((PadicInt.toZModPow n (M 1 0)).val) • Q := by
    rw [hPdef, galPoint_tateVec]
    exact addHom_mulVec_single n (tateVec W l e n) (tateVec_torsion W l e n) M 0
  have hσQ : galPoint W σ Q
      = ((PadicInt.toZModPow n (M 0 1)).val) • P
        + ((PadicInt.toZModPow n (M 1 1)).val) • Q := by
    rw [hQdef, galPoint_tateVec]
    exact addHom_mulVec_single n (tateVec W l e n) (tateVec_torsion W l e n) M 1
  have hmain := galois_cyclotomic W (l ^ n) hm σ P Q hPt hQt
    ((PadicInt.toZModPow n (M 0 0)).val) ((PadicInt.toZModPow n (M 0 1)).val)
    ((PadicInt.toZModPow n (M 1 0)).val) ((PadicInt.toZModPow n (M 1 1)).val)
    hσP hσQ hgen ζ hζ
  rw [galRep_coe, ← hM]
  exact eq_pow_of_mul_eq hm hζ (det_modEq n M) hmain

/-! ## ★出典の紐付け(`.src`) -/

def det_cyclotomic_of_gen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式が円分指標であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
