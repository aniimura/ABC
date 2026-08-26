import ABC3.Found.GaloisRep.WeilDet

/-!
# Galois (G5) 第 204 ブロック —— **★★★★★★★★★円分指標の段**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★`det ρ = 円分指標` の中身

基底 `P, Q` に `σ` が `σP = aP + cQ`・`σQ = bP + dQ` と効くとする。★このとき

    σ(e_m(P,Q)) = e_m(σP, σQ) = e_m(P,Q)^{ad − bc}

(第 193 の Galois 同変性と第 203 の行列式の公式)。★★`e_m(P,Q)` が `μ_m` を生成すれば、
任意の `ζ ∈ μ_m` は `e_m(P,Q)^k` と書けるので、`k` 乗して

    σ ζ = ζ^{ad − bc}

が出る。★★★これが**円分指標**である。

### ★引き算を避ける

指数を `ℕ` に保つため、`σ ζ · ζ^{bc} = ζ^{ad}` の形で述べる。

### ★★生成性(`hgen`)が非退化性である

「`ζ^m = 1` なら `ζ = e_m(P,Q)^k`」——これは `e_m(P,Q)` が原始 `m` 乗根であること、
すなわち**非退化性**である。★第 197 で `F(E)^{E[n]} = [n]^*F(E)` の 1 つに絞れており、
それは `deg[n] = n²`、すなわち EDS 恒等式(**mathlib 自身の TODO**)に帰着している。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `galois_cyclotomic` | ★★★★★★★★★**`σ ζ · ζ^{bc} = ζ^{ad}`** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]
  [IsAlgClosed L] [CharZero L] (W : WeierstrassCurve K)
  [((W.baseChange L).toAffine).IsElliptic]

/-- ★★★★★★★★★**円分指標の段**——`σ ζ = ζ^{ad−bc}`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 193 の Galois 同変性と第 203 の行列式の公式を、`hgen`(生成性 = 非退化性)で
`μ_m` 全体に広げる。★★指数を `ℕ` に保つため `σ ζ · ζ^{bc} = ζ^{ad}` の形。 -/
theorem galois_cyclotomic (m : ℕ) (hm : 1 ≤ m) (σ : L ≃ₐ[K] L)
    (P Q : (W.baseChange L).toAffine.Point)
    (hP : m • P = 0) (hQ : m • Q = 0) (a b c d : ℕ)
    (hσP : galPoint W σ P = a • P + c • Q)
    (hσQ : galPoint W σ Q = b • P + d • Q)
    (hgen : ∀ ζ : L, ζ ^ m = 1 →
      ∃ k : ℕ, (weilPairingVal (W.baseChange L).toAffine m P Q) ^ k = ζ)
    (ζ : L) (hζ : ζ ^ m = 1) :
    σ ζ * ζ ^ (b * c) = ζ ^ (a * d) := by
  have hbase : σ (weilPairingVal (W.baseChange L).toAffine m P Q)
      * (weilPairingVal (W.baseChange L).toAffine m P Q) ^ (b * c)
      = (weilPairingVal (W.baseChange L).toAffine m P Q) ^ (a * d) := by
    rw [weilPairingVal_galPoint W m hm σ P Q, hσP, hσQ]
    exact weilPairingVal_det (W.baseChange L).toAffine m hm P Q hP hQ a b c d
  obtain ⟨k, hk⟩ := hgen ζ hζ
  rw [← hk, map_pow, ← pow_mul, ← pow_mul]
  calc (σ (weilPairingVal (W.baseChange L).toAffine m P Q)) ^ k
        * (weilPairingVal (W.baseChange L).toAffine m P Q) ^ (k * (b * c))
      = ((σ (weilPairingVal (W.baseChange L).toAffine m P Q))
          * (weilPairingVal (W.baseChange L).toAffine m P Q) ^ (b * c)) ^ k := by
        rw [mul_pow, ← pow_mul, Nat.mul_comm k (b * c)]
    _ = ((weilPairingVal (W.baseChange L).toAffine m P Q) ^ (a * d)) ^ k := by rw [hbase]
    _ = (weilPairingVal (W.baseChange L).toAffine m P Q) ^ (k * (a * d)) := by
        rw [← pow_mul, Nat.mul_comm]

/-! ## ★出典の紐付け(`.src`) -/

def galois_cyclotomic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式が円分指標であること——円分指標の段)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
