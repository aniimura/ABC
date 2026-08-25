import ABC3.Found.GaloisRep.CollGroupK

/-!
# Galois (G6) 第 280 ブロック —— **★★★★★★★★領域の分離と原点近傍の単射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★付値を使わずに領域を分ける

群法則(第 279)には 3 点の `x` 座標の相異性が要る。原点近傍(`1 − p ∈ 𝔪`)の点は
`X` が極を持つので、他の領域の点とは自動的に異なる——これを**付値を持ち出さずに**示す:

    X·(1−p)² = XE、  XE は単元、  1 − p は非単元

★★★もし `X = algebraMap x` なら `XE = x(1−p)² ∈ 𝔪` となり、`XE` が単元であることに
反する。**「単元は真のイデアルに入らない」だけで済む**(`tateXK_not_mem_range`)。

★★`XE = p + (1−p)²(…)` で `p ≡ 1` は単元、`(1−p)²(…) ∈ 𝔪` だから `XE` は単元である。
第 272 の `ne_of_isUnit_of_mem`(単元 対 `I` の元)と同じ骨である。

## ★★★★★★★★原点近傍の中での単射性

`X` だけでは `±` を除いてしか決まらない。曲線の式の差から

    (Y − Y')(Y + Y' + X) = 0

第 1 の枝では `Y = Y'` なので `z = −X/Y`(第 276)も一致し、`tateZ_inj` で `u = u'`。
★★★**`z` が `X` と `Y` から復元できる**のが効いている——母数を商の形に選んだ理由。

## ★★相異性の全体像

| 組 | 相異性の根拠 |
|---|---|
| 単元 対 環帯 | `ne_of_isUnit_of_mem`(第 272) |
| 原点近傍 対 その他 | `tateXK_ne_of_origin`(**本ブロック**) |
| 単元どうし | `tateXpair_ne_of_units`(第 273) |
| 環帯どうし | `tateXpair_ne_of_ne`(第 260) |
| 原点近傍どうし | `tateXK_inj_origin`(**本ブロック**) |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isUnit_tateXpairE` | ★★★★原点近傍では `XE` は単元 |
| `tateXK_not_mem_range` | ★★★★★★★★**原点近傍の `X` は `R` の像に入らない** |
| `tateXK_ne_of_origin` | ★★★★★★★★**領域が違えば座標が違う** |
| `tate_equationK` | ★★★★★★`K` の水準の方程式(係数を明示) |
| `tateXK_inj_origin` | ★★★★★★★★**原点近傍での `X` の単射性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★原点近傍では `XE` は単元 -/

/-- ★★★★**原点近傍では `XE` は単元**——`XE ≡ p ≡ 1`。 -/
theorem isUnit_tateXpairE [IsAdicComplete I R] (p r q : R) (hq : q ∈ I) (hp : 1 - p ∈ I) :
    IsUnit (tateXpairE p r q hq) := by
  have hpU : IsUnit p := isUnit_of_one_sub_mem hp
  have hmem : (1 - p) ^ 2 * (tateXtail p q hq + (tateXterm r + tateXtail r q hq)
      - 2 * evalAdic (sigmaSeries 1) q hq) ∈ I := by
    refine Ideal.mul_mem_right _ _ ?_
    rw [pow_succ]
    exact Ideal.mul_mem_left _ _ hp
  rw [tateXpairE]
  exact isUnit_add_mem hpU hmem

/-! ## ★★★★★★★★領域の分離 -/

section Separate

variable {K : Type} [IsAdicComplete I R] [Nontrivial R] [Field K] [Algebra R K]

/-- ★★★★★★★★**原点近傍の `X` は `R` の像に入らない**——極を持つから。

★付値は要らない:`X·(1−p)² = XE` で `XE` は単元、`1−p` は非単元。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateXK_not_mem_range (hinj : Function.Injective (algebraMap R K))
    (p r q : R) (hq : q ∈ I) (hp : 1 - p ∈ I) (hp0 : algebraMap R K (1 - p) ≠ 0) (x : R) :
    (tateXK p r q hq : K) ≠ algebraMap R K x := by
  intro h
  rw [tateXK] at h
  have h2 : algebraMap R K (tateXpairE p r q hq) = algebraMap R K (x * (1 - p) ^ 2) := by
    rw [map_mul, map_pow]
    field_simp at h
    linear_combination h
  have h3 := hinj h2
  have h4 : tateXpairE p r q hq ∈ I := by
    rw [h3, pow_succ]
    exact Ideal.mul_mem_left _ _ (Ideal.mul_mem_left _ _ hp)
  exact not_isUnit_of_mem h4 (isUnit_tateXpairE p r q hq hp)

/-- ★★★★★★★★**領域が違えば座標が違う**(原点近傍 対 その他)。 -/
theorem tateXK_ne_of_origin (hinj : Function.Injective (algebraMap R K))
    (p r p' r' q : R) (hq : q ∈ I) (hp : 1 - p ∈ I) (hp0 : algebraMap R K (1 - p) ≠ 0)
    (hp' : IsUnit (1 - p')) :
    (tateXK p r q hq : K) ≠ tateXK p' r' q hq := by
  rw [tateXK_eq p' r' q hq hp']
  exact tateXK_not_mem_range hinj p r q hq hp hp0 _

end Separate

/-! ## ★★★★★★★★原点近傍の中での単射性 -/

section Inj

variable {K : Type} [IsAdicComplete I R] [Field K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★`K` の水準の方程式(係数を明示した形)。 -/
theorem tate_equationK (a w q : R) (hq : q ∈ I) (haw : a * w = q) (hw : IsUnit (1 - w))
    (ha : algebraMap R K (1 - a) ≠ 0) :
    tateYK (K := K) a w q hq ^ 2 + tateXK a w q hq * tateYK (K := K) a w q hq
      = tateXK a w q hq ^ 3 + algebraMap R K ((tateCurveAt q hq).a₄) * tateXK a w q hq
        + algebraMap R K ((tateCurveAt q hq).a₆) := by
  have h := congrArg (algebraMap R K) (tate_equationE a w q hq haw hw)
  simp only [map_add, map_mul, map_pow] at h
  rw [tateXK, tateYK]
  field_simp
  linear_combination (algebraMap R K (1 - a)) ^ 0 * h

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**原点近傍での `X` の単射性**(`±` を除く)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateXK_inj_origin (hinj : Function.Injective (algebraMap R K))
    (u u' q : R) (hq : q ∈ I) (hu : 1 - u ∈ I) (hu' : 1 - u' ∈ I)
    (hn : algebraMap R K (1 - u) ≠ 0) (hn' : algebraMap R K (1 - u') ≠ 0)
    (hX : tateXK u (wOf q u) q hq = (tateXK u' (wOf q u') q hq : K)) :
    u = u' ∨ tateYK (K := K) u' (wOf q u') q hq
      = -tateXK u (wOf q u) q hq - tateYK (K := K) u (wOf q u) q hq := by
  have e1 := tate_equationK (K := K) u (wOf q u) q hq
    (mul_wOf q u (isUnit_of_one_sub_mem hu)) (isUnit_one_sub (wOf_mem u hq)) hn
  have e2 := tate_equationK (K := K) u' (wOf q u') q hq
    (mul_wOf q u' (isUnit_of_one_sub_mem hu')) (isUnit_one_sub (wOf_mem u' hq)) hn'
  rw [← hX] at e2
  have hfac : (tateYK (K := K) u (wOf q u) q hq - tateYK (K := K) u' (wOf q u') q hq)
      * (tateYK (K := K) u (wOf q u) q hq + tateYK (K := K) u' (wOf q u') q hq
        + tateXK u (wOf q u) q hq) = 0 := by
    linear_combination e1 - e2
  rcases mul_eq_zero.1 hfac with h | h
  · refine Or.inl ?_
    have hY : tateYK (K := K) u (wOf q u) q hq = tateYK (K := K) u' (wOf q u') q hq :=
      sub_eq_zero.1 h
    have hz : algebraMap R K (tateZ u (wOf q u) q hq)
        = algebraMap R K (tateZ u' (wOf q u') q hq) := by
      rw [tateZ_map u (wOf q u) q hq hu hn, tateZ_map u' (wOf q u') q hq hu' hn', hX, hY]
    exact tateZ_inj u u' q hq hu hu' (hinj hz)
  · exact Or.inr (by linear_combination h)

end Inj

/-! ## ★出典の紐付け(`.src`) -/

def tateXK_not_mem_range.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍の X は R の像に入らない)",
    sectionId := "genell-def-3-3" }

def tateXK_inj_origin.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍での X の単射性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
