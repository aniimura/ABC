import ABC3.Found.GaloisRep.GalRepWitness

/-!
# スケルトン —— **行列式は円分指標(Weil 対)**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★(G5) に残った数学

(G3)(G4) は達成した——表現 `ρ_{E,l}` は `T_l E` への**本物の作用**の行列表示として
`Found/GaloisRep/GalRep.lean` に構成済みである。

★(G5) `FullImageData` に残っているのは `det_cyclotomic`:

    det ρ(σ) は 1 の `l^n` 乗根への σ の作用の指数である

★★これは **Weil 対**の内容である——`∧² T_l E ≅ ℤ_l(1)` そのもの。

## ★★★★証明が要求するもの

| 段 | 内容 | mathlib |
|---|---|---|
| (a) | `E[n]` 上の Weil 対 `e_n : E[n] × E[n] → μ_n` の構成 | ★**0 件** |
| (b) | 双線型性・交代性 | ★**0 件** |
| (c) | 非退化性(基底で原始 `n` 乗根になる) | ★**0 件** |
| (d) | Galois 同変性 `e_n(σP, σQ) = σ(e_n(P,Q))` | ★**0 件** |

★★(a) の古典的構成は「`div(f) = n(P) − n(O)` なる関数 `f`」と
「`g^n = f∘[n]`」を使い、**関数体上の因子**が要る。
★★★mathlib は `CoordinateRing`・`ClassGroup`・`toClass` を持つが、
**平行移動 `τ_Q` と乗法 `[n]` の関数体への引き戻し**が無い(2026-08-20 実測)。

★★★★見積もり: **50–150 ブロック**(因子の層から積むことになる)。

## ★★これも「壁」ではない

`Found/GaloisRep/` には既に `E[n] ≃+ (ℤ/n)²`(第 65–72)と `T_l E ≃+ ℤ_l²`(第 73–77)がある。
★Weil 対はその上に乗る**独立した層**である。
-/

namespace ABC3.Skeleton.GaloisRep

open ABC3.Meta ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve

/-! ## ★★★★★★行列式は円分指標 -/

/-- ★★★★★★**`det ρ(σ)` は円分指標である**(Weil 対の内容)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★表現 `galRep` は `Found/GaloisRep/GalRep.lean` に**構成済み**である。
★★残っているのはこの等式だけであり、`.needs` に (a)–(d) を書き下した。 -/
theorem det_galRep_eq_cyclotomic {K L : Type} [Field K] [DecidableEq K] [CharZero K]
    [Field L] [DecidableEq L] [Algebra K L] [IsAlgClosed L]
    (W : WeierstrassCurve K) (hell : W.IsElliptic) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L) (n : ℕ) (ζ : L) (hζ : ζ ^ (l ^ n) = 1) :
    σ ζ = ζ ^ ((PadicInt.toZModPow n
      ((galRep W l e σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).det)).val) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def det_galRep_eq_cyclotomic.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式が円分指標であること)",
    sectionId := "genell-thm-3-8" }

def det_galRep_eq_cyclotomic.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★★★**2026-08-20: 数学はすべて揃った**(第 178-204 ブロック)。★Weil 対 `e_n` の構成(第 178)、`e_n^n = 1`(第 179)、双線型性(第 184・191・195)、交代性(第 190・191)、Galois 同変性(第 192-194)、反対称性(第 195)、**行列式の公式 `e_n(aP+cQ, bP+dQ)·e_n(P,Q)^{bc} = e_n(P,Q)^{ad}`**(第 203)、**円分指標の段 `σ ζ · ζ^{bc} = ζ^{ad}`**(第 204)。★★当初の見積もり『50-150 ブロック、因子の層から積む』は**2 つの残件**に絞れた(0 ブロック)" 19,
    .implicitStep
      "★★★★★★★★**残件 (i): 非退化性 = `normEDS` が楕円列であること**。第 197 で非退化性は `F(E)^{E[n]} = [n]^*F(E)` の 1 つに絞れ、第 196 の Artin(`[F(E) : F(E)^{E[n]}] = n²`)と第 198 のモニック性(`Φ_n − c·ΨSq_n` は次数ちょうど `n²`)で挟むと、要るのは `x(nP) = Φ_n/ΨSq_n` だけになる。★これは分点多項式列が**楕円列**であることに帰着し、**mathlib 自身の TODO** である(`Mathlib/NumberTheory/EllipticDivisibilitySequence.lean` の `TODO: prove that normEDS satisfies IsEllDivSequence`、2026-08-20 実測)(30-80 ブロック、上流案件)" 19,
    .implicitStep
      "★★★★★★**残件 (ii): Lean 上の配線**。`galRep W l e σ` は `galTate` の基底 `e` での行列 `galMat` である(`Found/GaloisRep/GalRep.lean`)。★第 204 を当てるには、`T_l E` の基底 `e⁻¹(1,0)`・`e⁻¹(0,1)` の第 `n` 成分 `P, Q` が `E[l^n]` の基底になり、`σP = aP + cQ`・`σQ = bP + dQ` の `(a b; c d)` が `galMat mod l^n` に一致することが要る。★★要するに **`T_l E / l^n·T_l E ≅ E[l^n]`** と、`ℤ_l` の `E[l^n]` への作用が `PadicInt.toZModPow n` を経由すること。★★★射影 `T_l E → E[l^n]` は `Interface/GaloisRep/Torsion.lean` に定義済み、全射性に要る `E[l^{n+1}] → E[l^n]` の全射は**第 186 の `exists_nsmul_eq_point` で済んでいる**(5-15 ブロック)" 19,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, III.8(Weil 対の構成と性質)"
      (.absent "mathlib に Weil 対は 0 件(2026-08-20、`WeilPairing|weil_pairing` で全文検索して 0 件)") 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の Galois 同変性)" 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の双線型性)" 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対の非退化性)" 19,
    .otherPaper "GenEll" "Theorem 3.8(Weil 対が 1 の n 乗根に値を取ること)" 19,
    .implicitStep
      "★★★★2026-08-20 の展開: 本節点の下に部分グラフを立てた。葉は `Skeleton/GaloisRep/WeilFunctionField.lean` の 2 件——平行移動 τ_Q と乗法 [n] の関数体への引き戻しである" 19,
    .implicitStep
      "★★★★★`f_P`(div(f_P) = n(P) − n(O))は **Found に入った**(第 113 ブロック `xyIdeal_pow_isPrincipal`)。mathlib が群法則を類群経由で証明しているため `Point.toClass` と `ClassGroup.mk_eq_one_iff` がそのまま使えた。当初の「因子の層から積む」という見積もりは**下方修正**される(0 ブロック)" 19,
    .implicitStep
      "★★★mathlib は CoordinateRing・ClassGroup・toClass・分点多項式(Φ・ΨSq)を持つが、**平行移動 τ_Q と乗法 [n] の関数体への引き戻し**が無い(2026-08-20 実測)" 19,
    .implicitStep
      "★★★★E[n] ≃+ (ℤ/n)²(第 65-72)と T_l E ≃+ ℤ_l²(第 73-77)は Found/ に済んでいる。Weil 対はその上の独立した層である" 19 ]

end ABC3.Skeleton.GaloisRep
