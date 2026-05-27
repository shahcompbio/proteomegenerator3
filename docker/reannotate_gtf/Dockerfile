FROM golang:1.22-alpine AS builder

WORKDIR /build
COPY bin/reannotate_gtf.go .
RUN go build -ldflags="-s -w" -o reannotate_gtf reannotate_gtf.go

FROM alpine:3.20

LABEL org.opencontainers.image.source="https://github.com/kentsislab/proteomegenerator3"
LABEL org.opencontainers.image.description="reannotate_gtf: Re-annotate gffcompare GTFs with reference IDs"

COPY --from=builder /build/reannotate_gtf /usr/local/bin/reannotate_gtf

ENTRYPOINT ["reannotate_gtf"]
