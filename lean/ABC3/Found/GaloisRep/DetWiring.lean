import ABC3.Found.GaloisRep.TateWiring

/-!
# Galois (G5) 第 206 ブロック —— **★★★★★★行列成分の読み替え**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★配線の最後の 3 本

第 205 で `φ(α·v) = (toZModPow n α).val · φ(v)` が取れたので、
`galMat` の成分を `E[l^n]` の基底で読み替える段に進む。

### ★★★★★★1 本目——列を基底で読む

`M *ᵥ Pi.single j 1 = fun i => M i j` であり、これは
`(M 0 j)·single₀ + (M 1 j)·single₁` である。★したがって

    φ(M *ᵥ single j 1) = (toZModPow n (M 0 j)).val · φ(single₀)
                        + (toZModPow n (M 1 j)).val · φ(single₁)

★★`P := φ(single₀)`・`Q := φ(single₁)` と置けば、これが
`σP = aP + cQ`・`σQ = bP + dQ` そのものである。

### ★★★★★2 本目——行列式の合同

`det M = M₀₀M₁₁ − M₀₁M₁₀` に `toZModPow n`(環準同型)を当てて

    (toZModPow n det).val + B·C ≡ A·D    (mod l^n)

★`ZMod (l^n)` で計算してから `.val` に戻す。

### ★★★3 本目——指数は `l^n` を法として効く

`ζ^{l^n} = 1` なら `j ≡ k (mod l^n)` のとき `ζ^j = ζ^k`。
★これで第 204 の `σ ζ · ζ^{bc} = ζ^{ad}` を
`σ ζ = ζ^{(toZModPow n det).val}` に書き換えられる。

## ★★残るのは基底性だけ

`P, Q` が `E[l^n]` の**基底**であること——`T_l E ≅ ℤ_l²` の基底の第 `n` 成分が
`E[l^n] ≅ (ℤ/l^n)²` の基底になること。★射影の全射性に要る `E[l^{n+1}] → E[l^n]` は
**第 186 の `exists_nsmul_eq_point` で済んでいる**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `addHom_mulVec_single` | ★★★★★★**列を基底で読む** |
| `det_modEq` | ★★★★★**行列式の合同** |
| `pow_eq_pow_of_modEq` | ★★★指数は `l^n` を法として効く |
-/

namespace ABC3.Found.GaloisRep

universe u

/-- ★★★指数は `m` を法として効く。 -/
theorem pow_eq_pow_of_modEq {F : Type} [Field F] {m : ℕ} {ζ : F} (hζ : ζ ^ m = 1)
    {j k : ℕ} (h : j % m = k % m) : ζ ^ j = ζ ^ k := by
  conv_lhs => rw [← Nat.div_add_mod j m]
  conv_rhs => rw [← Nat.div_add_mod k m]
  rw [pow_add, pow_add, pow_mul, pow_mul, hζ, one_pow, one_pow, one_mul, one_mul, h]

/-- ★★★★★★**行列の列を `E[l^n]` の基底で読む**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`P := φ(single₀)`・`Q := φ(single₁)` と置けば `σP = aP + cQ` そのものである。 -/
theorem addHom_mulVec_single {A : Type u} [AddCommGroup A] {l : ℕ} [Fact l.Prime]
    (n : ℕ) (φ : (Fin 2 → ℤ_[l]) →+ A) (htor : ∀ v, (l ^ n) • φ v = 0)
    (M : Matrix (Fin 2) (Fin 2) ℤ_[l]) (j : Fin 2) :
    φ (Matrix.mulVec M (Pi.single j 1))
      = ((PadicInt.toZModPow n (M 0 j)).val) • φ (Pi.single 0 1)
        + ((PadicInt.toZModPow n (M 1 j)).val) • φ (Pi.single 1 1) := by
  have hcol : Matrix.mulVec M (Pi.single j 1)
      = (M 0 j) • (Pi.single 0 1 : Fin 2 → ℤ_[l])
        + (M 1 j) • (Pi.single 1 1 : Fin 2 → ℤ_[l]) := by
    funext i
    fin_cases i <;> simp [Matrix.mulVec, dotProduct, Pi.single_apply]
  rw [hcol, map_add, addHom_padic_smul_of_torsion n φ htor,
    addHom_padic_smul_of_torsion n φ htor]

/-- ★★★★★**行列式の合同**——`det + BC ≡ AD (mod l^n)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem det_modEq {l : ℕ} [Fact l.Prime] (n : ℕ) (M : Matrix (Fin 2) (Fin 2) ℤ_[l]) :
    ((PadicInt.toZModPow n M.det).val
        + (PadicInt.toZModPow n (M 0 1)).val * (PadicInt.toZModPow n (M 1 0)).val)
      % (l ^ n)
      = ((PadicInt.toZModPow n (M 0 0)).val * (PadicInt.toZModPow n (M 1 1)).val) % (l ^ n) := by
  have hne : NeZero (l ^ n) := ⟨by
    have h1 : 1 ≤ l ^ n := Nat.one_le_pow _ _ (Nat.pos_of_ne_zero (Fact.out : l.Prime).ne_zero)
    omega⟩
  refine (ZMod.natCast_eq_natCast_iff' _ _ _).1 ?_
  push_cast [ZMod.natCast_val, ZMod.cast_id]
  rw [Matrix.det_fin_two]
  simp only [map_sub, map_mul]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def addHom_mulVec_single.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——行列の列を基底で読むこと)",
    sectionId := "genell-thm-3-8" }

def det_modEq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——行列式の合同)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
