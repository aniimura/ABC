import ABC3.Found.GaloisRep.CyclotomicStep

/-!
# Galois (G5) 第 205 ブロック —— **★★★★★★Tate 加群から `E[l^n]` への配線**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★残件 (ii) の 2 本

第 204 の円分指標の段を `det_galRep_eq_cyclotomic` に繋ぐには、`T_l E` の基底の
第 `n` 成分が `E[l^n]` でどう動くかが要る。★必要なのは次の 2 本である。

### ★★★1 本目——射影は Galois 作用と可換

    tateProj n (galTate σ f) = galPoint σ (tateProj n f)

★`galTate` は成分ごとに `galPoint` を当てているので **`rfl`** で済む。

### ★★★★★★2 本目——`ℤ_l` の作用は `toZModPow` に落ちる

`φ := (E[l^n] への射影) ∘ e⁻¹ : (Fin 2 → ℤ_l) →+ E[l^n]` は**加法的**で、
値は `l^n` 捩れである。★このとき

    φ(α • v) = (toZModPow n α).val • φ(v)

★★証明は 3 行:`α − (toZModPow n α).val` は `toZModPow n` で消えるので
`l^n·γ`(第 78 の `padic_dvd_of_toZModPow`)、その分は `l^n` 捩れで消える。
★★★**`limTors` に `ℤ_l` 加群の構造を入れる必要は無い**——加法性と捩れだけで足りる。

## ★★これで残るのは行列成分の読み替えだけ

`P := φ(Pi.single 0 1)`・`Q := φ(Pi.single 1 1)` と置くと、
`galMat_apply`(`e (galTate σ x) = M *ᵥ e x`)から

    σP = (toZModPow n (M 0 0)).val • P + (toZModPow n (M 1 0)).val • Q

が出る(2 本目を `M *ᵥ Pi.single 0 1 = fun i => M i 0` に当てる)。
★★あとは `P, Q` が `E[l^n]` の基底であること(第 186 の全射性から)と、
`ad − bc = det M mod l^n` を合わせれば第 204 が当たる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateProj_galTate` | ★★★射影は Galois 作用と可換 |
| `addHom_padic_smul_of_torsion` | ★★★★★★**`ℤ_l` の作用は `toZModPow` に落ちる** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

universe u

/-- ★★★射影は Galois 作用と可換である。 -/
theorem tateProj_galTate {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L]
    [Algebra K L] (W : WeierstrassCurve K) (l n : ℕ) (σ : L ≃ₐ[K] L)
    (f : tateModule (W.baseChange L) l) :
    ((tateProj (W.baseChange L) l n (galTate W l σ f) : (W.baseChange L).toAffine.Point))
      = galPoint W σ ((tateProj (W.baseChange L) l n f : (W.baseChange L).toAffine.Point)) :=
  rfl

/-- ★★★★★★**加法写像が捩れに値を取るなら、`ℤ_l` の作用は `toZModPow` に落ちる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`limTors` に `ℤ_l` 加群の構造を入れる必要は無い——加法性と `l^n` 捩れだけで足りる。 -/
theorem addHom_padic_smul_of_torsion {A : Type u} [AddCommGroup A] {l : ℕ} [Fact l.Prime]
    (n : ℕ) (φ : (Fin 2 → ℤ_[l]) →+ A) (htor : ∀ v, (l ^ n) • φ v = 0)
    (α : ℤ_[l]) (v : Fin 2 → ℤ_[l]) :
    φ (α • v) = ((PadicInt.toZModPow n α).val) • φ v := by
  have hne : NeZero (l ^ n) := ⟨by
    have h1 : 1 ≤ l ^ n := Nat.one_le_pow _ _ (Nat.pos_of_ne_zero (Fact.out : l.Prime).ne_zero)
    omega⟩
  set a : ℕ := (PadicInt.toZModPow n α).val with ha
  have hz : PadicInt.toZModPow n (α - (a : ℤ_[l])) = 0 := by
    rw [map_sub, map_natCast, ha]
    simp
  obtain ⟨γ, hγ⟩ := padic_dvd_of_toZModPow l n _ hz
  have hsplit : α = (a : ℤ_[l]) + (l : ℤ_[l]) ^ n * γ := by
    rw [← hγ]; ring
  rw [hsplit, add_smul, map_add]
  have h1 : φ ((a : ℤ_[l]) • v) = a • φ v := by
    rw [Nat.cast_smul_eq_nsmul ℤ_[l] a v, map_nsmul]
  have h2 : φ (((l : ℤ_[l]) ^ n * γ) • v) = 0 := by
    have hs : ((l : ℤ_[l]) ^ n * γ) • v = (l ^ n) • (γ • v) := by
      rw [← Nat.cast_smul_eq_nsmul ℤ_[l] (l ^ n) (γ • v), Nat.cast_pow, mul_smul]
    rw [hs, map_nsmul, htor]
  rw [h1, h2, add_zero]

/-! ## ★出典の紐付け(`.src`) -/

def tateProj_galTate.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——射影が Galois 作用と可換であること)",
    sectionId := "genell-thm-3-8" }

def addHom_padic_smul_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——Z_l の作用が toZModPow に落ちること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
