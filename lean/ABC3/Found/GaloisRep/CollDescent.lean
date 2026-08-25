import ABC3.Found.GaloisRep.CollDivisibility

/-!
# Galois (G6) 第 253 ブロック —— **★★★★★★★3 変数の降下の骨格**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★2 変数のときとの唯一の違い

第 239 の降下は `ℤ[A,W] = ℤ[a][w]` で、`W^m ∣ P` を `X^m ∣ toPP P` に直し、
係数(**1 変数**多項式)が無限個の点で消えることから 0 を出した。

★3 変数では `ℤ[U,V,W] = ℤ[u][v][w]` となり、係数が **2 変数**になる。
そこで **「無限集合の直積の上で消える 2 変数多項式は 0」**(`poly2_eq_zero_of_infinite`)
を 1 枚はさむ。★これは 1 変数の議論を 2 回入れ子にするだけで、
まず `v` を動かして `ℂ[v]` の元として 0、次にその係数を `u` で動かして 0、とする。

## ★★★配線

| 写像 | 内容 |
|---|---|
| `toPP3 : ℤ[U,V,W] → ℤ[u][v][w]` | `U ↦ C(C X)`、`V ↦ C X`、`W ↦ X` |
| `ofPP3` | 逆向き |
| `evalUV u v : ℤ[U,V,W] → ℂ[w]` | `U ↦ C u`、`V ↦ C v`、`W ↦ X` |
| `eval2Z u v : ℤ[u][v] → ℂ` | 係数を `(u,v)` で評価 |

★`evalUV u v = (map (eval2Z u v)) ∘ toPP3` なので、
`evalUV u v P` の係数がすべて消えることと `toPP3 P` の係数が `(u,v)` で消えることが
同じになる(`coeff_evalUV_eq`)。

## ★★★★★★★到達点

    (∀ u v ∈ s, ∀ j < n+1, (evalUV u v P).coeff j = 0)  ⟹  W^{n+1} ∣ P

★あとは「`s` を上半平面の指数像に取り、係数の消滅を解析的な評価から出す」だけになる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `cW_pow_dvd_iff` | ★★★★`W^m ∣ P ↔ X^m ∣ toPP3 P` |
| `coeff_evalUV_eq` | ★★★★係数の対応 |
| `intPoly_eq_zero_of_infinite` | ★★★★1 変数の消滅 |
| `poly2_eq_zero_of_infinite` | ★★★★★★**2 変数の消滅** |
| `cW_pow_dvd_of_coeff` | ★★★★★★★**係数の消滅から整除性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real MvPolynomial

/-! ## ★★★★★`ℤ[U,V,W] ≅ ℤ[u][v][w]` -/

/-- ★`ℤ[u][v][w]`。 -/
abbrev PPP : Type := Polynomial (Polynomial (Polynomial ℤ))

/-- ★★★★★`ℤ[U,V,W] → ℤ[u][v][w]`(`W` を一番外側の変数に)。 -/
noncomputable def toPP3 : CollBase →+* PPP :=
  MvPolynomial.eval₂Hom (S₁ := PPP)
    (((Polynomial.C : Polynomial (Polynomial ℤ) →+* PPP).comp
      (Polynomial.C : Polynomial ℤ →+* Polynomial (Polynomial ℤ))).comp
      (Polynomial.C : ℤ →+* Polynomial ℤ))
    ![Polynomial.C (Polynomial.C Polynomial.X), Polynomial.C Polynomial.X, Polynomial.X]

/-- ★逆向き。 -/
noncomputable def ofPP3 : PPP →+* CollBase :=
  Polynomial.eval₂RingHom
    (Polynomial.eval₂RingHom
      (Polynomial.eval₂RingHom (MvPolynomial.C : ℤ →+* CollBase) cU) cV) cW

theorem toPP3_cU : toPP3 cU = Polynomial.C (Polynomial.C Polynomial.X) := by
  rw [toPP3, cU, MvPolynomial.eval₂Hom_X']; rfl

theorem toPP3_cV : toPP3 cV = Polynomial.C Polynomial.X := by
  rw [toPP3, cV, MvPolynomial.eval₂Hom_X']; rfl

theorem toPP3_cW : toPP3 cW = Polynomial.X := by
  rw [toPP3, cW, MvPolynomial.eval₂Hom_X']; rfl

theorem ofPP3_toPP3 : ofPP3.comp toPP3 = RingHom.id CollBase := by
  refine MvPolynomial.ringHom_ext ?_ ?_
  · intro r
    simp [toPP3, ofPP3]
  · intro i
    fin_cases i
    · show ofPP3 (toPP3 cU) = cU
      rw [toPP3_cU, ofPP3]; simp
    · show ofPP3 (toPP3 cV) = cV
      rw [toPP3_cV, ofPP3]; simp
    · show ofPP3 (toPP3 cW) = cW
      rw [toPP3_cW, ofPP3]; simp

/-- ★★★★**`W^m ∣ P` は `X^m ∣ toPP3 P` と同じ**。 -/
theorem cW_pow_dvd_iff (m : ℕ) (P : CollBase) :
    cW ^ m ∣ P ↔ Polynomial.X ^ m ∣ toPP3 P := by
  constructor
  · rintro ⟨c, hc⟩
    exact ⟨toPP3 c, by rw [hc, map_mul, map_pow, toPP3_cW]⟩
  · rintro ⟨g, hg⟩
    refine ⟨ofPP3 g, ?_⟩
    have h := congrArg ofPP3 hg
    rw [map_mul, map_pow] at h
    have hW : ofPP3 Polynomial.X = cW := by rw [ofPP3]; simp
    rw [hW] at h
    have hid : ofPP3 (toPP3 P) = P := congrArg (fun f => f P) ofPP3_toPP3
    rw [hid] at h
    exact h

/-! ## ★★★★分子を `w` の多項式と見る -/

/-- ★★`U ↦ u`・`V ↦ v`(定数)、`W ↦ X` として `ℤ[U,V,W]` を `ℂ[w]` に送る。 -/
noncomputable def evalUV (u v : ℂ) : CollBase →+* Polynomial ℂ :=
  MvPolynomial.eval₂Hom ((Polynomial.C : ℂ →+* Polynomial ℂ).comp (Int.castRingHom ℂ))
    ![Polynomial.C u, Polynomial.C v, Polynomial.X]

theorem evalUV_cU (u v : ℂ) : evalUV u v cU = Polynomial.C u := by
  rw [evalUV, cU, MvPolynomial.eval₂Hom_X']; rfl

theorem evalUV_cV (u v : ℂ) : evalUV u v cV = Polynomial.C v := by
  rw [evalUV, cV, MvPolynomial.eval₂Hom_X']; rfl

theorem evalUV_cW (u v : ℂ) : evalUV u v cW = Polynomial.X := by
  rw [evalUV, cW, MvPolynomial.eval₂Hom_X']; rfl

/-- ★★★★**`w` を代入すると `collEval u v w` に戻る**。 -/
theorem eval_evalUV (u v w : ℂ) (P : CollBase) :
    (evalUV u v P).eval w = collEval u v w P := by
  have hhom : (Polynomial.evalRingHom w).comp (evalUV u v) = collEval u v w := by
    refine MvPolynomial.ringHom_ext ?_ ?_
    · intro r
      simp [evalUV, collEval]
    · intro i
      fin_cases i
      · simp [evalUV, collEval, MvPolynomial.eval₂Hom_X']
      · simp [evalUV, collEval, MvPolynomial.eval₂Hom_X']
      · simp [evalUV, collEval, MvPolynomial.eval₂Hom_X']
  exact congrArg (fun f => f P) hhom

/-- ★`ℤ[u][v]` を `(u,v)` で評価する。 -/
noncomputable def eval2Z (u v : ℂ) : Polynomial (Polynomial ℤ) →+* ℂ :=
  Polynomial.eval₂RingHom (Polynomial.eval₂RingHom (Int.castRingHom ℂ) u) v

theorem evalUV_eq_map_toPP3 (u v : ℂ) :
    evalUV u v = (Polynomial.mapRingHom (eval2Z u v)).comp toPP3 := by
  refine MvPolynomial.ringHom_ext ?_ ?_
  · intro r
    simp [evalUV, toPP3, eval2Z]
  · intro i
    fin_cases i
    · show evalUV u v cU = _
      rw [evalUV_cU]
      show _ = (Polynomial.mapRingHom _) (toPP3 cU)
      rw [toPP3_cU]
      simp [eval2Z]
    · show evalUV u v cV = _
      rw [evalUV_cV]
      show _ = (Polynomial.mapRingHom _) (toPP3 cV)
      rw [toPP3_cV]
      simp [eval2Z]
    · show evalUV u v cW = _
      rw [evalUV_cW]
      show _ = (Polynomial.mapRingHom _) (toPP3 cW)
      rw [toPP3_cW]
      simp

/-- ★★★★**`evalUV u v` は `toPP3` の係数を `(u,v)` で評価したもの**。 -/
theorem coeff_evalUV_eq (u v : ℂ) (P : CollBase) (j : ℕ) :
    (evalUV u v P).coeff j = eval2Z u v ((toPP3 P).coeff j) := by
  have h : evalUV u v P = (Polynomial.mapRingHom (eval2Z u v)) (toPP3 P) :=
    congrArg (fun f : CollBase →+* Polynomial ℂ => f P) (evalUV_eq_map_toPP3 u v)
  rw [h]
  simp [Polynomial.coeff_map]

/-! ## ★★★★★★2 変数の消滅 -/

/-- ★★★★無限集合の上で消える `ℤ` 係数 1 変数多項式は 0。 -/
theorem intPoly_eq_zero_of_infinite (f : Polynomial ℤ) (s : Set ℂ) (hs : s.Infinite)
    (h : ∀ u ∈ s, Polynomial.eval₂ (Int.castRingHom ℂ) u f = 0) : f = 0 := by
  have hmap : f.map (Int.castRingHom ℂ) = 0 := by
    refine poly_eq_zero_of_infinite_zeros _ s hs ?_
    intro u hu
    rw [Polynomial.eval_map]
    exact h u hu
  have hcast : Function.Injective ((Int.castRingHom ℂ) : ℤ →+* ℂ) := by
    intro a b hab
    simpa using hab
  have hinj : Function.Injective (Polynomial.map (Int.castRingHom ℂ)) :=
    Polynomial.map_injective _ hcast
  exact hinj (by rw [hmap]; simp)

/-- ★★★★★★**無限集合の直積の上で消える 2 変数多項式は 0**。

★1 変数の議論を 2 回入れ子にするだけ——まず `v` を動かして `ℂ[v]` の元として 0、
次にその係数を `u` で動かして 0。 -/
theorem poly2_eq_zero_of_infinite (g : Polynomial (Polynomial ℤ)) (s : Set ℂ) (hs : s.Infinite)
    (h : ∀ u ∈ s, ∀ v ∈ s, eval2Z u v g = 0) : g = 0 := by
  refine Polynomial.ext fun j => ?_
  rw [Polynomial.coeff_zero]
  refine intPoly_eq_zero_of_infinite _ s hs ?_
  intro u hu
  have hgu : g.map (Polynomial.eval₂RingHom (Int.castRingHom ℂ) u) = 0 := by
    refine poly_eq_zero_of_infinite_zeros _ s hs ?_
    intro v hv
    rw [Polynomial.eval_map]
    exact h u hu v hv
  have hc := congrArg (fun p : Polynomial ℂ => p.coeff j) hgu
  simp only [Polynomial.coeff_map, Polynomial.coeff_zero] at hc
  exact hc

/-! ## ★★★★★★★係数の消滅から整除性 -/

/-- ★★★★★★★**`(u,v)` を動かして `W^{n+1} ∣ P` を出す**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem cW_pow_dvd_of_coeff (n : ℕ) (P : CollBase) (s : Set ℂ) (hs : s.Infinite)
    (h : ∀ u ∈ s, ∀ v ∈ s, ∀ j, j < n + 1 → (evalUV u v P).coeff j = 0) :
    cW ^ (n + 1) ∣ P := by
  rw [cW_pow_dvd_iff, Polynomial.X_pow_dvd_iff]
  intro j hj
  refine poly2_eq_zero_of_infinite _ s hs ?_
  intro u hu v hv
  rw [← coeff_evalUV_eq]
  exact h u hu v hv j hj

/-! ## ★出典の紐付け(`.src`) -/

def poly2_eq_zero_of_infinite.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——2 変数多項式の消滅)",
    sectionId := "genell-def-3-3" }

def cW_pow_dvd_of_coeff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——係数の消滅から W の整除性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
